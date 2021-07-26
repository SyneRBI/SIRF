/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2019 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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
\ingroup MR
\brief Implementation file for extended Gadgetron functionality classes.

\author Evgueni Ovtchinnikov
\author SyneRBI
*/
#include "sirf/iUtilities/DataHandle.h"
#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"
#include "sirf/Gadgetron/gadgetron_client.h"

using namespace gadgetron;
using namespace sirf;

GTConnector::GTConnector()
{
	sptr_con_ = gadgetron::shared_ptr < GadgetronClientConnector >
		(new GadgetronClientConnector);
}
GadgetronClientConnector& GTConnector::operator()()
{
	return *sptr_con_.get();
}
gadgetron::shared_ptr<GadgetronClientConnector> GTConnector::sptr()
{
	return sptr_con_;
}

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
	for (gh = gadgets_.begin(); gh != gadgets_.end(); ++gh) {
		if (sirf::iequals(gh->get()->id(), id))
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
	for (gh = readers_.begin(); gh != readers_.end(); ++gh)
		xml_script += gh->get()->gadget().xml() + '\n';
	for (gh = writers_.begin(); gh != writers_.end(); ++gh)
		xml_script += gh->get()->gadget().xml() + '\n';
	for (gh = gadgets_.begin(); gh != gadgets_.end(); ++gh) {
		const GadgetHandle* ptr_gh = gh->get();
		xml_script += ptr_gh->gadget().xml(ptr_gh->id()) + '\n';
//		xml_script += gh->get()->gadget().xml() + '\n';
	}
	if (endgadget_.get())
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
	uint32_t nacq = acquisitions.number();
	//std::cout << nacq << " acquisitions" << std::endl;
	sptr_acqs_ = acquisitions.new_acquisitions_container();
	if (nacq < 1)
		return;

	ISMRMRD::Acquisition acq_tmp;
	std::string config = xml();

	// quick fix: checking if AcquisitionFinishGadget is needed (= running old Gadgetron)
	shared_ptr<MRAcquisitionData> sptr_acqs = acquisitions.new_acquisitions_container();
	{
		GTConnector conn;
		conn().register_reader(GADGET_MESSAGE_ISMRMRD_ACQUISITION,
			shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientAcquisitionMessageCollector(sptr_acqs)));
		for (int nt = 0; nt < N_TRIALS; nt++) {
//			std::cout << "connection attempt " << nt << '\n';
			try {
				conn().connect(host_, port_);
				conn().send_gadgetron_configuration_script(config);
				conn().send_gadgetron_parameters(acquisitions.acquisitions_info());
				acquisitions.get_acquisition(0, acq_tmp);
				conn().send_ismrmrd_acquisition(acq_tmp);
				conn().send_gadgetron_close();
				conn().wait();
				break;
			}
			catch (...) {
				if (connection_failed(nt))
					THROW("Server running Gadgetron not accessible");
			}
		}
	}

	uint32_t na = sptr_acqs->number();
	//std::cout << na << " acquisitions processed\n";
	if (na < 1) {
		// old Gadgetron is running, have to append AcquisitionFinishGadget to the chain
		gadgetron::shared_ptr<AcquisitionFinishGadget>
			endgadget(new AcquisitionFinishGadget);
		set_endgadget(endgadget);
		config = xml();
	}

	GTConnector conn;
	conn().register_reader(GADGET_MESSAGE_ISMRMRD_ACQUISITION,
		shared_ptr<GadgetronClientMessageReader>
		(new GadgetronClientAcquisitionMessageCollector(sptr_acqs_)));
	conn().connect(host_, port_);
	conn().send_gadgetron_configuration_script(config);
	conn().send_gadgetron_parameters(acquisitions.acquisitions_info());
	for (uint32_t i = 0; i < nacq; i++) {
		acquisitions.get_acquisition(i, acq_tmp);
		conn().send_ismrmrd_acquisition(acq_tmp);
	}
	conn().send_gadgetron_close();
	conn().wait();
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
ImagesProcessor::process(const GadgetronImageData& images)
{
	std::string config = xml();
	GTConnector conn;
	sptr_images_ = images.new_images_container();
	if (dicom_)
		conn().register_reader(GADGET_MESSAGE_DICOM_WITHNAME,
			shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientBlobMessageReader(prefix_, "dcm")));
	else
		conn().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(sptr_images_)));
	for (int nt = 0; nt < N_TRIALS; nt++) {
		try {
			conn().connect(host_, port_);
			conn().send_gadgetron_configuration_script(config);
			for (unsigned int i = 0; i < images.number(); i++) {
				if (dicom_)
					conn().send_wrapped_image(*images.image_wrap(i).abs());
				else
					conn().send_wrapped_image(images.image_wrap(i));
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
			if (sirf::iequals(mc.as_str(attr, i), "image")) {
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
MRAcquisitionModel::set_up(gadgetron::shared_ptr<MRAcquisitionData> sptr_ac, 
			gadgetron::shared_ptr<GadgetronImageData> sptr_ic)
{
	if( sptr_ac->number() ==0 )
		throw LocalisedException("Please dont use an empty acquisition template.", __FILE__, __LINE__);

	if(sptr_ac->get_trajectory_type() == ISMRMRD::TrajectoryType::CARTESIAN)
		this->sptr_enc_ = std::make_shared<sirf::CartesianFourierEncoding>();
	else if(sptr_ac->get_trajectory_type() == ISMRMRD::TrajectoryType::OTHER)
	{
		ASSERT(sptr_ac->get_trajectory_dimensions()>0, "You should set a type ISMRMRD::TrajectoryType::OTHER trajectory before calling the calculate method with dimension > 0.");
	#ifdef GADGETRON_TOOLBOXES_AVAILABLE
	#warning "We compile the non-cartesian code in GADGETRON_X"
		this->sptr_enc_ = std::make_shared<sirf::RPEFourierEncoding>();
	#else
		throw std::runtime_error("Non-cartesian reconstruction is not supported, but your file contains ISMRMRD::TrajectoryType::OTHER data.");
	#endif
	}
	else
		throw std::runtime_error("Only cartesian or OTHER type of trajectory are available.");

	sptr_acqs_ = sptr_ac;
	set_image_template(sptr_ic);
}

void
MRAcquisitionModel::fwd(GadgetronImageData& ic, CoilSensitivitiesVector& cc,
	MRAcquisitionData& ac)
{
    GadgetronImagesVector images_channelresolved;
    cc.forward(images_channelresolved, ic);

    ac.sort();
    std::vector<KSpaceSubset> kspace_sorting = ac.get_kspace_sorting();

    if( kspace_sorting.size() != images_channelresolved.number() )
        throw LocalisedException("Number of images does not match number of acquisition data bins  ", __FILE__, __LINE__);

    for( unsigned int i=0; i<images_channelresolved.number(); ++i)
    {
        ImageWrap iw = images_channelresolved.image_wrap(i);
        CFImage* ptr_img = static_cast<CFImage*>(iw.ptr_image());

        auto tag_img = KSpaceSubset::get_tag_from_img(*ptr_img);

        sirf::AcquisitionsVector subset;
        KSpaceSubset::SetType idx_set;
        for(int j=0; j<kspace_sorting.size(); ++j)
        {
            if(tag_img == kspace_sorting[j].get_tag())
            {
                idx_set = kspace_sorting[j].get_idx_set();
                ac.get_subset(subset, idx_set);
                break;
            }
        }

        if(subset.number() == 0)
            throw LocalisedException("You didn't find rawdata corresponding to your image in the acquisition data.", __FILE__, __LINE__);

        this->sptr_enc_->forward(subset, *ptr_img);
        ac.set_subset(subset, idx_set); //assume forward does not reorder the acquisitions
    }
    ac.sort();
}

void 
MRAcquisitionModel::bwd(GadgetronImageData& ic, const CoilSensitivitiesVector& cc,
    const MRAcquisitionData& ac)
{
    GadgetronImagesVector iv;
    iv.set_meta_data(ac.acquisitions_info());

    auto sort_idx = ac.get_kspace_order();

    for(int i=0; i<sort_idx.size(); ++i)
    {
        sirf::AcquisitionsVector subset;
        ac.get_subset(subset, sort_idx[i]);
        
		CFImage* img_ptr = new CFImage();
		ImageWrap iw(ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT, img_ptr);// God I trust this!
		this->sptr_enc_->backward(*img_ptr, subset);

        iv.append(iw);
	}

    cc.backward(ic, iv);
    ic.set_up_geom_info();
}
