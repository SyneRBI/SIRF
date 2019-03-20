/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

/*!
\file
\ingroup xGadgetron Client
\brief Utilities for data exchange between SIRF and Gadgetron server.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef GADGETRON_CLIENT
#define GADGETRON_CLIENT

#include <chrono>
#include <condition_variable>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <thread>

#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>

#include <ismrmrd/dataset.h>
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>

#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"

enum GadgetronMessageID {
	GADGET_MESSAGE_INT_ID_MIN = 0,
	GADGET_MESSAGE_CONFIG_FILE = 1,
	GADGET_MESSAGE_CONFIG_SCRIPT = 2,
	GADGET_MESSAGE_PARAMETER_SCRIPT = 3,
	GADGET_MESSAGE_CLOSE = 4,
	GADGET_MESSAGE_TEXT = 5,
	GADGET_MESSAGE_INT_ID_MAX = 999,
	GADGET_MESSAGE_EXT_ID_MIN = 1000,
	GADGET_MESSAGE_ISMRMRD_ACQUISITION = 1008,
	GADGET_MESSAGE_CLOUD_JOB = 1013,
	GADGET_MESSAGE_GADGETCLOUD_JOB = 1014,
	GADGET_MESSAGE_DICOM_WITHNAME = 1018,
	GADGET_MESSAGE_DEPENDENCY_QUERY = 1019,
	GADGET_MESSAGE_ISMRMRD_IMAGE = 1022,
	GADGET_MESSAGE_RECONDATA = 1023,
	GADGET_MESSAGE_EXT_ID_MAX = 4096
};

namespace sirf {

	struct GadgetMessageIdentifier {
		uint16_t id;
	};

	struct GadgetMessageConfigurationFile {
		char configuration_file[1024];
	};

	struct GadgetMessageScript {
		uint32_t script_length;
	};

	class GadgetronClientException : public std::exception {
	public:
		GadgetronClientException(std::string msg) : msg_(msg) {}
		virtual ~GadgetronClientException() throw() {}
		virtual const char* what() const throw()
		{
			return msg_.c_str();
		}
	protected:
		std::string msg_;
	};

	/**
	\brief Abstract base class for receiving messages from Gadgetron server.
	*/
	class GadgetronClientMessageReader {
	public:
		virtual ~GadgetronClientMessageReader() {}
		/**
		Function must be implemented to read a specific message.
		*/
		virtual void read(boost::asio::ip::tcp::socket* s) = 0;
	};

	/**
	\brief Class for accumulating acquisitions sent by Gadgetron server.
	*/
	class GadgetronClientAcquisitionMessageCollector :
		public GadgetronClientMessageReader {
	public:
		GadgetronClientAcquisitionMessageCollector
			(gadgetron::shared_ptr<MRAcquisitionData> ptr_acqs) : ptr_acqs_(ptr_acqs) {}
		virtual ~GadgetronClientAcquisitionMessageCollector() {}

		virtual void read(boost::asio::ip::tcp::socket* stream);

	private:
		gadgetron::shared_ptr<MRAcquisitionData> ptr_acqs_;
	};

	/**
	\brief Class for accumulating images sent by Gadgetron server.
	*/
	class GadgetronClientImageMessageCollector :
		public GadgetronClientMessageReader {
	public:
		GadgetronClientImageMessageCollector
			(gadgetron::shared_ptr<GadgetronImageData> ptr_images) : ptr_images_(ptr_images) {}
		virtual ~GadgetronClientImageMessageCollector() {}

		template <typename T>
		void read_data_attributes
			(ISMRMRD::Image<T>* ptr, const ISMRMRD::ImageHeader& h, void** ptr_ptr,
			boost::asio::ip::tcp::socket* stream)
		{
			ISMRMRD::Image < T >* ptr_im = new ISMRMRD::Image < T > ;
			*ptr_ptr = (void*)ptr_im;
			ISMRMRD::Image<T>& im = *ptr_im;
			im.setHead(h);
			//im.setImageType(ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE);

			//Read meta attributes
			typedef unsigned long long size_t_type;
			size_t_type meta_attrib_length;
			boost::asio::read
				(*stream, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));
			if (meta_attrib_length > 0) {
				std::string meta_attrib(meta_attrib_length, 0);
				boost::asio::read(*stream, boost::asio::buffer(const_cast < char* >
					(meta_attrib.c_str()), meta_attrib_length));
				im.setAttributeString(meta_attrib);
			}
			//Read image data
			boost::asio::read
				(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
		}

		virtual void read(boost::asio::ip::tcp::socket* stream);

	private:
		gadgetron::shared_ptr<GadgetronImageData> ptr_images_;
	};

	/**
	\brief Class for communicating with Gadgetron server.
	*/
	class GadgetronClientConnector {
	public:
		GadgetronClientConnector() : socket_(0), timeout_ms_(2000)
		{}
		virtual ~GadgetronClientConnector()
		{
			if (socket_) {
				socket_->close();
				delete socket_;
			}
		}

		void set_timeout(unsigned int t)
		{
			timeout_ms_ = t;
		}

		void read_task();

		void wait()
		{
			reader_thread_.join();
		}

		void connect(std::string hostname, std::string port);

		void send_gadgetron_close();

		void send_gadgetron_configuration_file(std::string config_xml_name);

		void send_gadgetron_configuration_script(std::string xml_string);

		void  send_gadgetron_parameters(std::string xml_string);

		void send_ismrmrd_acquisition(ISMRMRD::Acquisition& acq);

		template<typename T>
		void send_ismrmrd_image(ISMRMRD::Image<T>* ptr_im)
		{
			ISMRMRD::Image<T>& im = *ptr_im;
			if (!socket_)
				throw GadgetronClientException("Invalid socket.");

			GadgetMessageIdentifier id;
			id.id = GADGET_MESSAGE_ISMRMRD_IMAGE;

			boost::asio::write
				(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
			boost::asio::write
				(*socket_,
				boost::asio::buffer(&im.getHead(), sizeof(ISMRMRD::ImageHeader)));

			size_t meta_attrib_length = im.getAttributeStringLength();
			std::string meta_attrib(meta_attrib_length + 1, 0);
			im.getAttributeString(meta_attrib);

			//std::cout << "attributes:" << std::endl << meta_attrib << std::endl;

			if (meta_attrib_length > 0) {
				size_t l = meta_attrib.find("</ismrmrdMeta>") + std::strlen("</ismrmrdMeta>");
				//std::cout << meta_attrib_length << ' ' << l << '\n';
				meta_attrib.erase(l);
			}

			boost::asio::write
				(*socket_, boost::asio::buffer(&meta_attrib_length, sizeof(size_t)));
			boost::asio::write
				(*socket_, boost::asio::buffer(meta_attrib.c_str(), meta_attrib_length));

			boost::asio::write
				(*socket_, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
		}

		void send_wrapped_image(ImageWrap& iw)
		{
			IMAGE_PROCESSING_SWITCH(iw.type(), send_ismrmrd_image, iw.ptr_image());
		}

		void register_reader
			(unsigned short slot, gadgetron::shared_ptr<GadgetronClientMessageReader> r)
		{
			readers_[slot] = r;
		}

	protected:
		typedef
			std::map < unsigned short, gadgetron::shared_ptr<GadgetronClientMessageReader> >
			maptype;

		GadgetronClientMessageReader* find_reader(unsigned short r);

		boost::asio::io_service io_service;
		boost::asio::ip::tcp::socket* socket_;
		boost::thread reader_thread_;
		maptype readers_;
		unsigned int timeout_ms_;
	};

}

#endif
