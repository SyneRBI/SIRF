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
\ingroup Gadgetron Extensions
\brief Implementation file for extended Gadgetron functionality classes.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#include "sirf/iUtilities/DataHandle.h"
#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/gadgetron_x.h"

using namespace gadgetron;
using namespace sirf;

static bool
connection_failed(int nt)
{
	std::cout << "connection to gadgetron server failed";
	if (nt < N_TRIALS - 1) {
		std::cout << ", trying again...\n";
		return false;
	}
	else {
		std::cout << std::endl;
		return true;
	}
}

static void
check_gadgetron_connection(std::string host, std::string port)
{
	ImagesProcessor ip;
	ip.set_host(host);
	ip.set_port(port);
	//std::cout << "checking connection...\n";
	for (int nt = 0; nt < N_TRIALS; nt++) {
		try {
			ip.check_connection();
			//std::cout << "ok\n";
			break;
		}
		catch (...) {
			if (connection_failed(nt))
				THROW("Connection to Gadgetron server lost, check Gadgetron output");
		}
	}
}

shared_ptr<aGadget> 
GadgetChain::gadget_sptr(std::string id)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
	std::list<shared_ptr<GadgetHandle> >::iterator gh;
#else
	typename std::list<shared_ptr<GadgetHandle> >::iterator gh;
#endif
	for (gh = gadgets_.begin(); gh != gadgets_.end(); gh++) {
		if (boost::iequals(gh->get()->id(), id))
			return gh->get()->gadget_sptr();
	}
	return shared_ptr<aGadget>();
}

std::string 
GadgetChain::xml() const 
{
	std::string xml_script("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	xml_script += "<gadgetronStreamConfiguration xsi:schemaLocation=";
	xml_script += "\"http://gadgetron.sf.net/gadgetron gadgetron.xsd\"\n";
	xml_script += "xmlns=\"http://gadgetron.sf.net/gadgetron\"\n";
	xml_script += "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n\n";

#if defined(_MSC_VER) && _MSC_VER < 1900
	std::list<shared_ptr<GadgetHandle> >::const_iterator gh;
#else
	typename std::list<shared_ptr<GadgetHandle> >::const_iterator gh;
#endif
	for (gh = readers_.begin(); gh != readers_.end(); gh++)
		xml_script += gh->get()->gadget().xml() + '\n';
	for (gh = writers_.begin(); gh != writers_.end(); gh++)
		xml_script += gh->get()->gadget().xml() + '\n';
	for (gh = gadgets_.begin(); gh != gadgets_.end(); gh++) {
		const GadgetHandle* ptr_gh = gh->get();
		xml_script += ptr_gh->gadget().xml(ptr_gh->id()) + '\n';
//		xml_script += gh->get()->gadget().xml() + '\n';
	}
	xml_script += endgadget_->xml() + '\n';
	xml_script += "</gadgetronStreamConfiguration>\n";

	return xml_script;
}

/*
The next three methods:

	AcquisitionsProcessor::process
	ImagesProcessor::process
	ImagesReconstructor::process

contain code snippets from
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

void 
AcquisitionsProcessor::process(MRAcquisitionData& acquisitions) 
{
	std::string config = xml();
	GTConnector conn;
	uint32_t nacq = 0;
	nacq = acquisitions.number();
	//std::cout << nacq << " acquisitions" << std::endl;
	ISMRMRD::Acquisition acq_tmp;
	sptr_acqs_ = acquisitions.new_acquisitions_container();

	conn().register_reader(GADGET_MESSAGE_ISMRMRD_ACQUISITION,
		shared_ptr<GadgetronClientMessageReader>
		(new GadgetronClientAcquisitionMessageCollector(sptr_acqs_)));
	for (int nt = 0; nt < N_TRIALS; nt++) {
		try {
			conn().connect(host_, port_);
			conn().send_gadgetron_configuration_script(config);
			conn().send_gadgetron_parameters(acquisitions.acquisitions_info());
			for (uint32_t i = 0; i < nacq; i++) {
				acquisitions.get_acquisition(i, acq_tmp);
				conn().send_ismrmrd_acquisition(acq_tmp);
			}
			conn().send_gadgetron_close();
			conn().wait();
			break;
		}
		catch (...) {
			if (connection_failed(nt))
				THROW("Server running Gadgetron not accessible");
		}
	}
	check_gadgetron_connection(host_, port_);
}

void 
ImagesReconstructor::process(MRAcquisitionData& acquisitions)
{
	//check_gadgetron_connection(host_, port_);

	std::string config = xml();
	//std::cout << "config:\n" << config << std::endl;

	uint32_t nacquisitions = 0;
	nacquisitions = acquisitions.number();
	//std::cout << nacquisitions << " acquisitions" << std::endl;
	ISMRMRD::Acquisition acq_tmp;

	GTConnector conn;
	//std::cout << "connecting to port " << port_ << "...\n";
	sptr_images_.reset(new GadgetronImagesVector);
	sptr_images_->set_meta_data(acquisitions.acquisitions_info());
	conn().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
		shared_ptr<GadgetronClientMessageReader>
		(new GadgetronClientImageMessageCollector(sptr_images_)));
	for (int nt = 0; nt < N_TRIALS; nt++) {
		try {
			conn().connect(host_, port_);
			conn().send_gadgetron_configuration_script(config);
			conn().send_gadgetron_parameters(acquisitions.acquisitions_info());
			for (uint32_t i = 0; i < nacquisitions; i++) {
				acquisitions.get_acquisition(i, acq_tmp);
				conn().send_ismrmrd_acquisition(acq_tmp);
			}
			conn().send_gadgetron_close();
			conn().wait();
			break;
		}
		catch (...) {
			if (connection_failed(nt))
				THROW("Server running Gadgetron not accessible");
		}
	}
	check_gadgetron_connection(host_, port_);
	sptr_images_->sort();
    // Add meta data to the image
    sptr_images_->set_meta_data(acquisitions.acquisitions_info());
}

void 
ImagesProcessor::process(GadgetronImageData& images)
{
	std::string config = xml();
	GTConnector conn;
	sptr_images_ = images.new_images_container();
	conn().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
		shared_ptr<GadgetronClientMessageReader>
		(new GadgetronClientImageMessageCollector(sptr_images_)));
	for (int nt = 0; nt < N_TRIALS; nt++) {
		try {
			conn().connect(host_, port_);
			conn().send_gadgetron_configuration_script(config);
			for (unsigned int i = 0; i < images.number(); i++) {
				ImageWrap& iw = images.image_wrap(i);
				conn().send_wrapped_image(iw);
			}
			conn().send_gadgetron_close();
			conn().wait();
			break;
		}
		catch (...) {
			if (connection_failed(nt))
				THROW("Server running Gadgetron not accessible");
		}
	}
	check_gadgetron_connection(host_, port_);
}

void
ImagesProcessor::check_connection()
{
	std::string config = xml();
	GTConnector conn;
	shared_ptr<GadgetronImageData> sptr_images(new GadgetronImagesVector);
	GadgetronImageData& images = *sptr_images_;
	conn().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
		shared_ptr<GadgetronClientMessageReader>
		(new GadgetronClientImageMessageCollector(sptr_images)));
	conn().connect(host_, port_);
	conn().send_gadgetron_configuration_script(config);
	ISMRMRD::Image<float>* ptr_im = new ISMRMRD::Image<float>(128, 128, 1);
	memset(ptr_im->getDataPtr(), 0, ptr_im->getDataSize());
	ImageWrap iw(ptr_im->getDataType(), ptr_im);
	conn().send_wrapped_image(iw);
	conn().send_gadgetron_close();
	conn().wait();
}

void MRAcquisitionModel::check_data_role(const GadgetronImageData& ic)
{
	for (int j = 0; j < ic.number(); j++) {
		const ImageWrap& iw = ic.image_wrap(j);
		std::string atts = iw.attributes();
		int atts_size = atts.size();
		if (atts_size < 1)
			continue;
		ISMRMRD::MetaContainer mc;
		ISMRMRD::deserialize(atts.c_str(), mc);
		const char* attr = "GADGETRON_DataRole";
		size_t l = mc.length(attr);
		bool ok = false;
		std::string value;
		for (int i = 0; i < l; i++) {
			if (boost::iequals(mc.as_str(attr, i), "image")) {
				ok = true;
				break;
			}
			if (i)
				value += " ";
			value += mc.as_str(attr, i);
		}
		if (!ok) {
			std::string msg("MRAcquisitionModel cannot use ");
			msg += "image data with GADGETRON_DataRole = ";
			msg += value;
			throw LocalisedException(msg.c_str(), __FILE__, __LINE__);
		}
	}
}

void
MRAcquisitionModel::fwd(GadgetronImageData& ic, CoilSensitivitiesContainer& cc, 
	MRAcquisitionData& ac)
{
	if (cc.items() < 1)
		throw LocalisedException
		("coil sensitivity maps not found", __FILE__, __LINE__);
	for (unsigned int i = 0, a = 0; i < ic.number(); i++) {
		ImageWrap& iw = ic.image_wrap(i);
		CoilData& csm = cc(i%cc.items());
		fwd(iw, csm, ac, a);
	}
}

void 
MRAcquisitionModel::bwd(GadgetronImageData& ic, CoilSensitivitiesContainer& cc, 
	MRAcquisitionData& ac)
{
	ic.set_meta_data(ac.acquisitions_info());
	if (cc.items() < 1)
		throw LocalisedException
		("coil sensitivity maps not found", __FILE__, __LINE__);
	for (unsigned int i = 0, a = 0; a < ac.number(); i++) {
		CoilData& csm = cc(i%cc.items());
		ImageWrap iw(sptr_imgs_->image_wrap(i));
		bwd(iw, csm, ac, a);
		ic.append(iw);
	}
}

/*
The next two methods:

	MRAcquisitionModel::fwd_
	MRAcquisitionModel::bwd_

contain code snippets from ISMRMRD/utilities/generate_cartesian_shepp_logan.cpp
by Michael S. Hansen

ISMRMRD SOFTWARE LICENSE JULY 2013

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

template< typename T>
void 
MRAcquisitionModel::fwd_(ISMRMRD::Image<T>* ptr_img, CoilData& csm,
	MRAcquisitionData& ac, unsigned int& off)
{
	ISMRMRD::Image<T>& img = *ptr_img;

	std::string par;
	ISMRMRD::IsmrmrdHeader header;
	//par = ac.acquisitions_info();
	par = sptr_acqs_->acquisitions_info();
	ISMRMRD::deserialize(par.c_str(), header);
	ISMRMRD::Encoding e = header.encoding[0];
	ISMRMRD::Acquisition acq; // (acq_);
	//sptr_acqs_->get_acquisition(0, acq);
	for (unsigned int i = 0; i < sptr_acqs_->number(); i++) {
		sptr_acqs_->get_acquisition(i, acq);
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
			break;
	}

	//int readout = e.encodedSpace.matrixSize.x;
	unsigned int nx = e.reconSpace.matrixSize.x;
	//unsigned int ny = e.reconSpace.matrixSize.y;
	//unsigned int nz = e.reconSpace.matrixSize.z;
	unsigned int ny = e.encodedSpace.matrixSize.y;
	unsigned int nz = e.encodedSpace.matrixSize.z;
	unsigned int nc = acq.active_channels();
	unsigned int readout = acq.number_of_samples();

	std::vector<size_t> dims;
	dims.push_back(readout);
	dims.push_back(ny);
	dims.push_back(nz);
	dims.push_back(nc);

	ISMRMRD::NDArray<complex_float_t> ci(dims);
	memset(ci.getDataPtr(), 0, ci.getDataSize());

	for (unsigned int c = 0; c < nc; c++) {
		for (unsigned int z = 0; z < nz; z++) {
			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++) {
					uint16_t xout = x + (readout - nx) / 2;
					complex_float_t zi = (complex_float_t)img(x, y, z);
					complex_float_t zc = csm(x, y, z, c);
					ci(xout, y, z, c) = zi * zc;
				}
			}
		}
	}

	memset((void*)acq.getDataPtr(), 0, acq.getDataSize());

	fft3c(ci);

	int y = 0;
	for (;;){
		sptr_acqs_->get_acquisition(off + y, acq);
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
			break;
		y++;
	}
	for (;;) {
		sptr_acqs_->get_acquisition(off + y, acq);
		int yy = acq.idx().kspace_encode_step_1;
		int zz = acq.idx().kspace_encode_step_2;
		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int s = 0; s < readout; s++) {
				acq.data(s, c) = ci(s, yy, zz, c);
			}
		}
		ac.append_acquisition(acq);
		y++;
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
			break;
	}
	off += y;

}

template< typename T>
void 
MRAcquisitionModel::bwd_(ISMRMRD::Image<T>* ptr_im, CoilData& csm,
	MRAcquisitionData& ac, unsigned int& off)
{
	ISMRMRD::Image<T>& im = *ptr_im;

	std::string par;
	ISMRMRD::IsmrmrdHeader header;
	par = ac.acquisitions_info();
	ISMRMRD::deserialize(par.c_str(), header);
	ISMRMRD::Encoding e = header.encoding[0];
	ISMRMRD::Acquisition acq;
	//sptr_acqs_->get_acquisition(0, acq);
	for (unsigned int i = 0; i < ac.number(); i++) {
		ac.get_acquisition(i, acq);
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
			break;
	}

	unsigned int nx = e.reconSpace.matrixSize.x;
	//unsigned int ny = e.reconSpace.matrixSize.y;
	//unsigned int nz = e.reconSpace.matrixSize.z;
	unsigned int ny = e.encodedSpace.matrixSize.y;
	unsigned int nz = e.encodedSpace.matrixSize.z;
	unsigned int nc = acq.active_channels();
	unsigned int readout = acq.number_of_samples();

	std::vector<size_t> dims;
	dims.push_back(readout);
	dims.push_back(ny);
	dims.push_back(nz);
	dims.push_back(nc);

	ISMRMRD::NDArray<complex_float_t> ci(dims);
	memset(ci.getDataPtr(), 0, ci.getDataSize());
	int y = 0;
	for (;;){
		ac.get_acquisition(off + y, acq);
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
			break;
		y++;
	}
	for (;;) {
		ac.get_acquisition(off + y, acq);
		int yy = acq.idx().kspace_encode_step_1;
		int zz = acq.idx().kspace_encode_step_2;
		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int s = 0; s < readout; s++) {
				ci(s, yy, zz, c) = acq.data(s, c);
			}
		}
		y++;
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
			break;
	}
	off += y;
	ifft3c(ci);

	T* ptr = im.getDataPtr();
	T s;
	memset(ptr, 0, im.getDataSize());
	long long int i = 0;
	for (unsigned int c = 0; c < nc; c++) {
		i = 0;
		for (unsigned int z = 0; z < nz; z++) {
			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++, i++) {
					uint16_t xout = x + (readout - nx) / 2;
					complex_float_t zi = ci(xout, y, z, c);
					complex_float_t zc = csm(x, y, z, c);
					xGadgetronUtilities::convert_complex(std::conj(zc) * zi, s);
					ptr[i] += s;
				}
			}
		}
	}

}

