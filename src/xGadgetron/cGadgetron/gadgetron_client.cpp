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

This file contains code snippets from
Gadgetron/apps/clients/gadgetron_ismrmrd_client/gadgetron_ismrmrd_client.cpp
by Michael S. Hansen

GADGETRON SOFTWARE LICENSE V1.0, NOVEMBER 2011

PERMISSION IS HEREBY GRANTED, FREE OF CHARGE, TO ANY PERSON OBTAINING
A COPY OF THIS SOFTWARE AND ASSOCIATED DOCUMENTATION FILES (THE
"SOFTWARE"), TO DEAL IN THE SOFTWARE WITHOUT RESTRICTION, INCLUDING
WITHOUT LIMITATION THE RIGHTS TO USE, COPY, MODIFY, MERGE, PUBLISH,
DISTRIBUTE, SUBLICENSE, AND/OR SELL COPIES OF THE SOFTWARE, AND TO
PERMIT PERSONS TO WHOM THE SOFTWARE IS FURNISHED TO DO SO, SUBJECT TO
THE FOLLOWING CONDITIONS:

THE ABOVE COPYRIGHT NOTICE, THIS PERMISSION NOTICE, AND THE LIMITATION
OF LIABILITY BELOW SHALL BE INCLUDED IN ALL COPIES OR REDISTRIBUTIONS
OF SUBSTANTIAL PORTIONS OF THE SOFTWARE.

SOFTWARE IS BEING DEVELOPED IN PART AT THE NATIONAL HEART, LUNG, AND BLOOD
INSTITUTE, NATIONAL INSTITUTES OF HEALTH BY AN EMPLOYEE OF THE FEDERAL
GOVERNMENT IN THE COURSE OF HIS OFFICIAL DUTIES. PURSUANT TO TITLE 17,
SECTION 105 OF THE UNITED STATES CODE, THIS SOFTWARE IS NOT SUBJECT TO
COPYRIGHT PROTECTION AND IS IN THE PUBLIC DOMAIN. EXCEPT AS CONTAINED IN
THIS NOTICE, THE NAME OF THE AUTHORS, THE NATIONAL HEART, LUNG, AND BLOOD
INSTITUTE (NHLBI), OR THE NATIONAL INSTITUTES OF HEALTH (NIH) MAY NOT
BE USED TO ENDORSE OR PROMOTE PRODUCTS DERIVED FROM THIS SOFTWARE WITHOUT
SPECIFIC PRIOR WRITTEN PERMISSION FROM THE NHLBI OR THE NIH.THE SOFTWARE IS
PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*!
\file
\ingroup SIRF Gadgetron client
\brief Implementation file for SIRF Gadgetron client.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#include <boost/asio.hpp>

#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/gadgetron_client.h"

using namespace gadgetron;
using namespace sirf;

void
GadgetronClientAcquisitionMessageCollector::read(boost::asio::ip::tcp::socket* stream)
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

void 
GadgetronClientImageMessageCollector::read(boost::asio::ip::tcp::socket* stream)
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
	}
	else {
		throw GadgetronClientException("Invalid image data type");
	}
}

void 
GadgetronClientConnector::read_task()
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

void 
GadgetronClientConnector::connect(std::string hostname, std::string port)
{
	boost::asio::ip::tcp::resolver resolver(io_service);
	boost::asio::ip::tcp::resolver::query 
		query(boost::asio::ip::tcp::v4(), hostname.c_str(), port.c_str());
	boost::asio::ip::tcp::resolver::iterator 
		endpoint_iterator = resolver.resolve(query);
	boost::asio::ip::tcp::resolver::iterator end;

	socket_ = new boost::asio::ip::tcp::socket(io_service);
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

	reader_thread_ = std::thread([&](){this->read_task();});
}

void 
GadgetronClientConnector::send_gadgetron_close()
{
	if (!socket_) {
		throw GadgetronClientException("Invalid socket.");
	}
	GadgetMessageIdentifier id;
	id.id = GADGET_MESSAGE_CLOSE;
	boost::asio::write
		(*socket_, boost::asio::buffer(&id, sizeof(GadgetMessageIdentifier)));
}

void 
GadgetronClientConnector::send_gadgetron_configuration_file(std::string config_xml_name)
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

void 
GadgetronClientConnector::send_gadgetron_configuration_script(std::string xml_string)
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

void 
GadgetronClientConnector::send_gadgetron_parameters(std::string xml_string)
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

void 
GadgetronClientConnector::send_ismrmrd_acquisition(ISMRMRD::Acquisition& acq)
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

GadgetronClientMessageReader* 
GadgetronClientConnector::find_reader(unsigned short r)
{
	GadgetronClientMessageReader* ret = 0;
	maptype::iterator it = readers_.find(r);
	if (it != readers_.end()) {
		ret = it->second.get();
	}
	return ret;
}

