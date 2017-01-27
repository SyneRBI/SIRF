#ifndef GADGETRON_CLIENT
#define GADGETRON_CLIENT

#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>

#include <fstream>
#include <streambuf>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <exception>
#include <map>
#include <thread>
#include <chrono>
#include <condition_variable>

using boost::asio::ip::tcp;

#include "gadgetron_data_containers.h"

enum GadgetronMessageID {
	GADGET_MESSAGE_INT_ID_MIN = 0,
	GADGET_MESSAGE_CONFIG_FILE = 1,
	GADGET_MESSAGE_CONFIG_SCRIPT = 2,
	GADGET_MESSAGE_PARAMETER_SCRIPT = 3,
	GADGET_MESSAGE_CLOSE = 4,
	GADGET_MESSAGE_TEXT = 5,
	GADGET_MESSAGE_INT_ID_MAX = 999,
	GADGET_MESSAGE_EXT_ID_MIN = 1000,
	GADGET_MESSAGE_ACQUISITION = 1008, /**< DEPRECATED */
	GADGET_MESSAGE_NEW_MEASUREMENT = 1002, /**< DEPRECATED */
	GADGET_MESSAGE_END_OF_SCAN = 1003, /**< DEPRECATED */
	GADGET_MESSAGE_IMAGE_CPLX_FLOAT = 1004, /**< DEPRECATED */
	GADGET_MESSAGE_IMAGE_REAL_FLOAT = 1005, /**< DEPRECATED */
	GADGET_MESSAGE_IMAGE_REAL_USHORT = 1006, /**< DEPRECATED */
	GADGET_MESSAGE_EMPTY = 1007, /**< DEPRECATED */
	GADGET_MESSAGE_ISMRMRD_ACQUISITION = 1008,
	GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT = 1009, /**< DEPRECATED */
	GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT = 1010, /**< DEPRECATED */
	GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT = 1011, /**< DEPRECATED */
	GADGET_MESSAGE_DICOM = 1012, /**< DEPRECATED */
	GADGET_MESSAGE_CLOUD_JOB = 1013,
	GADGET_MESSAGE_GADGETCLOUD_JOB = 1014,
	GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT = 1015, /**< DEPRECATED */
	GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT = 1016, /**< DEPRECATED */
	GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT = 1017, /**< DEPRECATED */
	GADGET_MESSAGE_DICOM_WITHNAME = 1018,
	GADGET_MESSAGE_DEPENDENCY_QUERY = 1019,
	GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_SHORT = 1020, /**< DEPRECATED */
	GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_SHORT = 1021, /**< DEPRECATED */
	GADGET_MESSAGE_ISMRMRD_IMAGE = 1022,
	GADGET_MESSAGE_RECONDATA = 1023,
	GADGET_MESSAGE_EXT_ID_MAX = 4096
};

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

class GadgetronClientMessageReader {
public:
	virtual ~GadgetronClientMessageReader() {}
	/**
	Function must be implemented to read a specific message.
	*/
	virtual void read(tcp::socket* s) = 0;
};

class GadgetronClientAcquisitionMessageCollector : 
	public GadgetronClientMessageReader {
public:
	GadgetronClientAcquisitionMessageCollector
		(boost::shared_ptr<AcquisitionsContainer> ptr_acqs) : ptr_acqs_(ptr_acqs) {}
	virtual ~GadgetronClientAcquisitionMessageCollector() {}

	virtual void read(tcp::socket* stream)
	{
		ISMRMRD::Acquisition acq;
		ISMRMRD::AcquisitionHeader h;
		boost::asio::read
			(*stream, boost::asio::buffer(&h, sizeof(ISMRMRD::AcquisitionHeader)));
		acq.setHead(h);
		unsigned long trajectory_elements =
			acq.getHead().trajectory_dimensions * acq.getHead().number_of_samples;
		unsigned long data_elements =
			acq.getHead().active_channels * acq.getHead().number_of_samples;

		if (trajectory_elements) {
			boost::asio::read
				(*stream, boost::asio::buffer
				(&acq.getTrajPtr()[0], sizeof(float)*trajectory_elements));
		}
		if (data_elements) {
			boost::asio::read
				(*stream,
				boost::asio::buffer
				(&acq.getDataPtr()[0], 2 * sizeof(float)*data_elements));
		}

		ptr_acqs_->append_acquisition(acq);
	}

private:
	boost::shared_ptr<AcquisitionsContainer> ptr_acqs_;
};

class GadgetronClientImageMessageCollector : 
	public GadgetronClientMessageReader {
public:
	GadgetronClientImageMessageCollector
		(boost::shared_ptr<ImagesContainer> ptr_images) : ptr_images_(ptr_images) {}
	virtual ~GadgetronClientImageMessageCollector() {}

	template <typename T>
	void read_data_attributes
		(ISMRMRD::Image<T>* ptr, const ISMRMRD::ImageHeader& h, void** ptr_ptr, 
		tcp::socket* stream)
	{
		ISMRMRD::Image < T >* ptr_im = new ISMRMRD::Image < T > ;
		*ptr_ptr = (void*)ptr_im;
		ISMRMRD::Image<T>& im = *ptr_im;
		im.setHead(h);
		im.setImageType(ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE);

		//Read meta attributes
		typedef unsigned long long size_t_type;
		size_t_type meta_attrib_length;
		boost::asio::read
			(*stream, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));
		if (meta_attrib_length>0) {
			std::string meta_attrib(meta_attrib_length, 0);
			boost::asio::read(*stream, boost::asio::buffer(const_cast<char*>
				(meta_attrib.c_str()), meta_attrib_length));
			im.setAttributeString(meta_attrib);
		}
		//Read image data
		boost::asio::read
			(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
	}

	virtual void read(tcp::socket* stream)
	{
		//Read the image headerfrom the socket
		ISMRMRD::ImageHeader h;
		boost::asio::read
			(*stream, boost::asio::buffer(&h, sizeof(ISMRMRD::ImageHeader)));
		void* ptr = 0;
		IMAGE_PROCESSING_SWITCH
			(h.data_type, read_data_attributes, ptr, h, &ptr, stream);
		if (ptr) {
			ptr_images_->append(h.data_type, ptr);
			ptr_images_->count(h.image_index);
		}
		else {
			throw GadgetronClientException("Invalid image data type");
		}
	}

private:
	boost::shared_ptr<ImagesContainer> ptr_images_;
};

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

	void read_task()
	{
		if (!socket_) {
			throw GadgetronClientException("Unable to create socket.");
		}

		GadgetMessageIdentifier id;
		while (socket_->is_open()) {
			try {
				boost::asio::read
					(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));

				if (id.id == GADGET_MESSAGE_CLOSE) {
					break;
				}

				GadgetronClientMessageReader* r = find_reader(id.id);

				if (!r) {
					std::cout << "Message received with ID: " << id.id << std::endl;
					throw GadgetronClientException("Unknown Message ID");
				}
				else {
					r->read(socket_);
				}
			}
			catch (...) {
				std::cout << "Input stream has terminated" << std::endl;
				return;
			}
		}
	}

	void wait() 
	{
		reader_thread_.join();
	}

	void connect(std::string hostname, std::string port)
	{
		tcp::resolver resolver(io_service);
		tcp::resolver::query query(tcp::v4(), hostname.c_str(), port.c_str());
		tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
		tcp::resolver::iterator end;

		socket_ = new tcp::socket(io_service);
		if (!socket_) {
			throw GadgetronClientException("Unable to create socket.");
		}

		std::condition_variable cv;
		std::mutex cv_m;

		boost::system::error_code error = boost::asio::error::host_not_found;
		std::thread t([&](){
			//TODO:
			//For newer versions of Boost, we should use
			//   boost::asio::connect(*socket_, iterator);
			while (error && endpoint_iterator != end) {
				socket_->close();
				socket_->connect(*endpoint_iterator++, error);
			}
			cv.notify_all();
		});

		{
			std::unique_lock<std::mutex> lk(cv_m);
			if (std::cv_status::timeout == 
				cv.wait_until(lk, std::chrono::system_clock::now() + 
				std::chrono::milliseconds(timeout_ms_))) {
				socket_->close();
			}
		}

		t.join();

		if (error)
			throw GadgetronClientException("Error connecting using socket.");

		reader_thread_ = 
			boost::thread(boost::bind(&GadgetronClientConnector::read_task, this));
	}

	void send_gadgetron_close() 
	{
		if (!socket_) {
			throw GadgetronClientException("Invalid socket.");
		}
		GadgetMessageIdentifier id;
		id.id = GADGET_MESSAGE_CLOSE;
		boost::asio::write
			(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
	}

	void send_gadgetron_configuration_file(std::string config_xml_name) 
	{
		if (!socket_) {
			throw GadgetronClientException("Invalid socket.");
		}

		GadgetMessageIdentifier id;
		id.id = GADGET_MESSAGE_CONFIG_FILE;

		GadgetMessageConfigurationFile ini;
		memset(&ini, 0, sizeof(GadgetMessageConfigurationFile));
		strncpy
			(ini.configuration_file, config_xml_name.c_str(), config_xml_name.size());

		boost::asio::write
			(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
		boost::asio::write
			(*socket_, 
			boost::asio::buffer(&ini, sizeof(GadgetMessageConfigurationFile)));
	}

	void send_gadgetron_configuration_script(std::string xml_string)
	{
		if (!socket_) {
			throw GadgetronClientException("Invalid socket.");
		}

		GadgetMessageIdentifier id;
		id.id = GADGET_MESSAGE_CONFIG_SCRIPT;

		GadgetMessageScript conf;
		conf.script_length = (uint32_t)xml_string.size() + 1;

		boost::asio::write
			(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
		boost::asio::write
			(*socket_, boost::asio::buffer(&conf, sizeof(GadgetMessageScript)));
		boost::asio::write
			(*socket_, boost::asio::buffer(xml_string.c_str(), conf.script_length));
	}

	void  send_gadgetron_parameters(std::string xml_string)
	{
		if (!socket_) {
			throw GadgetronClientException("Invalid socket.");
		}

		GadgetMessageIdentifier id;
		id.id = GADGET_MESSAGE_PARAMETER_SCRIPT;

		GadgetMessageScript conf;
		conf.script_length = (uint32_t)xml_string.size() + 1;

		boost::asio::write
			(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
		boost::asio::write
			(*socket_, boost::asio::buffer(&conf, sizeof(GadgetMessageScript)));
		boost::asio::write
			(*socket_, boost::asio::buffer(xml_string.c_str(), conf.script_length));
	}

	void send_ismrmrd_acquisition(ISMRMRD::Acquisition& acq)
	{
		if (!socket_)
			throw GadgetronClientException("Invalid socket.");

		GadgetMessageIdentifier id;
		id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

		boost::asio::write
			(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
		boost::asio::write
			(*socket_, 
			boost::asio::buffer(&acq.getHead(), sizeof(ISMRMRD::AcquisitionHeader)));

		unsigned long trajectory_elements = 
			acq.getHead().trajectory_dimensions*acq.getHead().number_of_samples;
		unsigned long data_elements = 
			acq.getHead().active_channels*acq.getHead().number_of_samples;

		if (trajectory_elements) {
			boost::asio::write
				(*socket_, boost::asio::buffer
				(&acq.getTrajPtr()[0], sizeof(float)*trajectory_elements));
		}
		if (data_elements) {
			boost::asio::write
				(*socket_, 
				boost::asio::buffer
				(&acq.getDataPtr()[0], 2 * sizeof(float)*data_elements));
		}
	}

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
		(unsigned short slot, boost::shared_ptr<GadgetronClientMessageReader> r) 
	{
		readers_[slot] = r;
	}

protected:
	typedef 
		std::map<unsigned short, boost::shared_ptr<GadgetronClientMessageReader> > 
		maptype;

	GadgetronClientMessageReader* find_reader(unsigned short r)
	{
		GadgetronClientMessageReader* ret = 0;
		maptype::iterator it = readers_.find(r);
		if (it != readers_.end()) {
			ret = it->second.get();
		}
		return ret;
	}

	boost::asio::io_service io_service;
	tcp::socket* socket_;
	boost::thread reader_thread_;
	maptype readers_;
	unsigned int timeout_ms_;
};

#endif
