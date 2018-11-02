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

#include <boost/asio.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

using boost::asio::ip::tcp;

#include "cgadgetron_shared_ptr.h"
#include "data_handle.h"
#include "gadgetron_x.h"

#include "encoding.h"

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
#ifdef _MSC_VER
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

#ifdef _MSC_VER
	std::list<shared_ptr<GadgetHandle> >::const_iterator gh;
#else
	typename std::list<shared_ptr<GadgetHandle> >::const_iterator gh;
#endif
	for (gh = readers_.begin(); gh != readers_.end(); gh++)
		xml_script += gh->get()->gadget().xml() + '\n';
	for (gh = writers_.begin(); gh != writers_.end(); gh++)
		xml_script += gh->get()->gadget().xml() + '\n';
	for (gh = gadgets_.begin(); gh != gadgets_.end(); gh++)
		xml_script += gh->get()->gadget().xml() + '\n';
	xml_script += endgadget_->xml() + '\n';
	xml_script += "</gadgetronStreamConfiguration>\n";

	return xml_script;
}

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
	sptr_images_.reset(new ImagesVector);
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
	sptr_images_->order();
}

void 
ImagesProcessor::process(MRImageData& images)
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
	shared_ptr<MRImageData> sptr_images(new ImagesVector);
	MRImageData& images = *sptr_images_;
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

void
MRAcquisitionModel::fwd(MRImageData& ic, CoilSensitivitiesContainer& cc, 
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
MRAcquisitionModel::bwd(MRImageData& ic, CoilSensitivitiesContainer& cc, 
	MRAcquisitionData& ac)
{
	if (cc.items() < 1)
		throw LocalisedException
		("coil sensitivity maps not found", __FILE__, __LINE__);
	ImageWrap iw(sptr_imgs_->image_wrap(0));
	for (unsigned int i = 0, a = 0; a < ac.number(); i++) {
		CoilData& csm = cc(i%cc.items());
		bwd(iw, csm, ac, a);
		ic.append(iw);
	}
}
/*
template< typename T>
void 
MRAcquisitionModel::fwd_(ISMRMRD::Image<T>* ptr_img, CoilData& csm,
	MRAcquisitionData& ac, unsigned int& off)
{
	ISMRMRD::Image<T>& img = *ptr_img;

	uint16_t const size_x = img.getMatrixSizeX(); 
	uint16_t const size_y = img.getMatrixSizeY(); 
	uint16_t const size_z = img.getMatrixSizeZ(); 
	uint16_t const size_dyn = img.getNumberOfChannels();	

	std::cout << " sx " << size_x << std::endl;
	std::cout << " sy " << size_y << std::endl;
	std::cout << " sz " << size_z << std::endl;
	std::cout << " sc " << size_dyn << std::endl;

	std::string par;
	ISMRMRD::IsmrmrdHeader header;
	par = sptr_acqs_->acquisitions_info();
	ISMRMRD::deserialize(par.c_str(), header);
	ISMRMRD::Encoding e = header.encoding[0];
		
	//int readout = e.encodedSpace.matrixSize.x;
	unsigned int nx = e.reconSpace.matrixSize.x;
	unsigned int ny = e.reconSpace.matrixSize.y;
	unsigned int nz = e.reconSpace.matrixSize.z;
	unsigned int nc = acq.active_channels();
	unsigned int readout = acq.number_of_samples();

	std::vector<size_t> dims;
	dims.push_back(readout);
	dims.push_back(ny);
	dims.push_back(nc);

	ISMRMRD::NDArray<complex_float_t> ci(dims);
	memset(ci.getDataPtr(), 0, ci.getDataSize());

	for (unsigned int c = 0; c < nc; c++) {
		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				uint16_t xout = x + (readout - nx) / 2;
				complex_float_t zi = (complex_float_t)img(x, y);
				complex_float_t zc = csm(x, y, 0, c);
				ci(xout, y, c) = zi * zc;
			}
		}
	}

	memset((void*)acq.getDataPtr(), 0, acq.getDataSize());

	fft2c(ci);

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
		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int s = 0; s < readout; s++) {
				acq.data(s, c) = ci(s, yy, c);
			}
		}
		ac.append_acquisition(acq);
		y++;
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
			break;
	}
	off += y;

}
*/
template< typename T>
void 
MRAcquisitionModel::fwd_(ISMRMRD::Image<T>* ptr_img, CoilData& csm,
	MRAcquisitionData& ac, unsigned int& off)
{
	ISMRMRD::Image<T>& img = *ptr_img;
	
	std::string par;
	ISMRMRD::IsmrmrdHeader header;

	par = sptr_acqs_->acquisitions_info();
	ISMRMRD::deserialize(par.c_str(), header);


	ISMRMRD::Encoding e = header.encoding[0];

	int readout_size = e.encodedSpace.matrixSize.x;
	
	unsigned int nx =img.getMatrixSizeX();//e.reconSpace.matrixSize.x;
	unsigned int ny =img.getMatrixSizeY();//e.reconSpace.matrixSize.y;
	unsigned int nz =img.getMatrixSizeZ();//e.reconSpace.matrixSize.z;



	if(  2 * img.getMatrixSizeX() != readout_size ) 
		throw LocalisedException("The image dimensions passed to the simulation are not half the size of the encoded space in the header. Readout os 2 is assumed.", __FILE__, __LINE__);

	if ( this->sptr_traj_->get_traj_type() == "Cartesian" )
	{
		if( img.getMatrixSizeY() != ny  )
			throw LocalisedException("Acquisition info contains ny not matching image dimension y", __FILE__, __LINE__);
		if( img.getMatrixSizeZ() != nz  )
			throw LocalisedException("Acquisition info contains nz not matching image dimension z", __FILE__, __LINE__);
	}
	// unsigned int nc = acq.active_channels();
	int coilmap_dims[4]; 
	csm.get_dim(coilmap_dims);
	unsigned int nc = coilmap_dims[3];


	std::vector<size_t> dims;
	dims.push_back(readout_size);
	dims.push_back(ny);
	dims.push_back(nz);
	dims.push_back(nc);
	
	ISMRMRD::NDArray<complex_float_t> ci(dims);
	memset(ci.getDataPtr(), 0, ci.getDataSize());

	// #pragma omp parallel for
	for (unsigned int c = 0; c < nc; c++) {
		for( unsigned int z = 0; z < nz; z++) {
			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++) {

					uint16_t xout = x + (readout_size - nx) / 2;
					complex_float_t zi = (complex_float_t)img(x, y, z);
					complex_float_t zc = csm(x, y, z, c);
					ci(xout, y, z, c) = zi * zc;
				}
			}
		}
	}



	ISMRMRD::NDArray< complex_float_t > k_data;
	
	std::string trajectory_type = this->sptr_traj_->get_traj_type();

	if( trajectory_type == "RPE" )
	{
		std::cout << "RPE Acquisition Process" << std::endl;

		RadialPhaseEncodingFFT RPE_FFT;
		auto traj = this->sptr_traj_->get_trajectory();
		RPE_FFT.set_trajectory( traj );
		RPE_FFT.SampleFourierSpace( ci );
		std::cout << "sampling done" << std::endl;

		k_data = RPE_FFT.get_k_data();

	}
	else if( trajectory_type == "" || trajectory_type == "Cartesian" ) 
	{
		std::cout << "Cartesian Acquisition Process" << std::endl;
		FullySampledCartesianFFT CartFFT;
		CartFFT.SampleFourierSpace( ci );
		
		k_data = CartFFT.get_k_data();
	}
	
	unsigned int const num_acq = sptr_acqs_->items(); 


	for( unsigned int i_acq = 0; i_acq < num_acq; i_acq++)
	{
		ISMRMRD::Acquisition acq; 
		sptr_acqs_->get_acquisition(i_acq, acq);

		unsigned int num_sampled_readout_pts = acq.number_of_samples();

		if( num_sampled_readout_pts < acq.number_of_samples() )			
			throw LocalisedException("The number of samples you try to acquire is larger than the volume dimension.", __FILE__, __LINE__);
			

		acq.resize(num_sampled_readout_pts, nc);
		memset((void*)acq.getDataPtr(), 0, acq.getDataSize());

		uint16_t const enc_step_1 = acq.getHead().idx.kspace_encode_step_1;
		uint16_t const enc_step_2 = acq.getHead().idx.kspace_encode_step_2;
		
		size_t const is_reverse = acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE)? 1: 0;

		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int s = 0; s < num_sampled_readout_pts; s++) {
				size_t const readout_access = is_reverse * (readout_size - 1 - s) + (1-is_reverse)*s;
				// acq.data(s, c) = k_data(s, enc_step_2, enc_step_1, c);
				acq.data(s, c) = k_data(readout_access, enc_step_2, enc_step_1, c);
			}
		}

		this->sptr_traj_->set_acquisition_trajectory(acq);
		acq.idx().contrast = img.getContrast();
	
		ac.append_acquisition(acq);

	}
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
	sptr_acqs_->get_acquisition(0, acq);

	unsigned int nx = e.reconSpace.matrixSize.x;
	unsigned int ny = e.reconSpace.matrixSize.y;
	unsigned int nc = acq.active_channels();
	unsigned int readout = acq.number_of_samples();

	std::vector<size_t> dims;
	dims.push_back(readout);
	dims.push_back(ny);
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
		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int s = 0; s < readout; s++) {
				ci(s, yy, c) = acq.data(s, c);
			}
		}
		y++;
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
			break;
	}
	off += y;
	ifft2c(ci);

	T* ptr = im.getDataPtr();
	T s;
	memset(ptr, 0, im.getDataSize());
	long long int i = 0;
	for (unsigned int c = 0; c < nc; c++) {
		i = 0;
		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++, i++) {
				uint16_t xout = x + (readout - nx) / 2;
				complex_float_t z = ci(xout, y, c);
				complex_float_t zc = csm(x, y, 0, c);
				xGadgetronUtilities::convert_complex(std::conj(zc) * z, s);
				ptr[i] += s;
			}
		}
	}

}

