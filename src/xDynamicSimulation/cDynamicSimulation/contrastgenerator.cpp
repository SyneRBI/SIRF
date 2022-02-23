/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "sirf/cDynamicSimulation/contrastgenerator.h"
#include "sirf/Reg/NiftyResample.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <omp.h>

//#include "Testing/auxiliary_testing_functions.h"

using namespace std;
using namespace stir;
using namespace sirf;
using namespace ISMRMRD;



AbstractContrastGenerator::AbstractContrastGenerator(const LabelVolume& tissue_labels, const std::string& filename_tissue_parameter_xml)
{
	this->tlm_ = TissueLabelMapper( tissue_labels, filename_tissue_parameter_xml );
	tlm_.map_labels_to_tissue_from_xml();

}

void AbstractContrastGenerator::replace_petmr_tissue_parameters(LabelType label, TissueParameter tiss_param)	
{
	this->tlm_.replace_petmr_tissue_parameters(label, tiss_param);
}

TissueParameter AbstractContrastGenerator::get_petmr_tissue_parameter( LabelType label ) 
{
	TissueParameterList all_tissues = this->tlm_.get_tissue_parameter_list();
	std::cout << "Number of different tissues: " << all_tissues.size()<< std::endl; 
	for( size_t i=0; i<all_tissues.size(); i++)
	{
		if( all_tissues[i].label_ == label)
			return all_tissues[i];
	}
	throw std::runtime_error("The label you asked for is apparently not given in the segmentation/XML combination you passed to the simulation.");
}		


MRContrastGenerator::MRContrastGenerator (const LabelVolume& tissue_labels, const std::string& filename_tissue_parameter_xml) :
AbstractContrastGenerator(tissue_labels, filename_tissue_parameter_xml)
{
}

void MRContrastGenerator::set_template_rawdata(const MRAcquisitionData& ad){

	this->sptr_acqu_ = std::move(ad.clone());
	this->hdr_ = this->sptr_acqu_->acquisitions_info().get_IsmrmrdHeader();
}

void MRContrastGenerator::set_rawdata_header(const ISMRMRD::IsmrmrdHeader& hdr)
{
	this->hdr_ = hdr;
}

sirf::GadgetronImagesVector& MRContrastGenerator::get_contrast_filled_volumes(bool const resample_output)
{
	if( resample_output == true )
		this->resample_to_template_image();

	return this->contrast_filled_volumes_;

}

void MRContrastGenerator::resample_to_template_image( void )
{
}

void MRContrastGenerator::map_contrast()
{
	this->tlm_.assign_tissues_to_labels();
	this->contrast_filled_volumes_.empty();
	if(sptr_acqu_ == nullptr)
		throw std::runtime_error("Your contrast template acquisition is not set.");

	this->contrast_filled_volumes_ = GadgetronImagesVector(*sptr_acqu_);
	
	std::vector < complex_float_t >	(*contrast_map_function)(std::shared_ptr<TissueParameter> const ptr_to_tiss_par, const ISMRMRD::IsmrmrdHeader& ismrmrd_hdr);

	ISMRMRD::SequenceParameters const sequ_par = *(this->hdr_.sequenceParameters);
	std::string const sequ_name = *(sequ_par.sequence_type);

	if(sequ_name.compare("Flash") == 0)
	{
		contrast_map_function = &map_flash_contrast;
	}
	else if (sequ_name.compare("Bssfp") == 0)
	{
		contrast_map_function = &map_bssfp_contrast;
	}
	else
	{
		std::stringstream error_msg_stream;
		error_msg_stream << "The header you read in requires a contrast which has not been implemented yet. ";
		error_msg_stream << "The demanded sequence type is: " << sequ_name << ". ";
		error_msg_stream << "Please give another rawdata header or write the contrast map yourself and add an else if to the map_contrast method.";
		throw std::runtime_error( error_msg_stream.str() );
	}

	TissueVector tissue_params = this->tlm_.get_segmentation_tissues();
	size_t const num_voxels = tissue_params.size();	

	std::vector<std::vector< complex_float_t> > contrast_vector;
	contrast_vector.resize(num_voxels);
	
	// #pragma omp parallel
	for (size_t i= 0; i<num_voxels; i++)
		contrast_vector[i] = contrast_map_function(tissue_params[i], this->hdr_);

	size_t const num_contrasts = contrast_vector[0].size();

	const int* segmentation_dims = this->tlm_.get_segmentation_dimensions();

	std::vector<size_t> data_size;
	
	for( int i_dim=0; i_dim< 8; i_dim++)
	{
		data_size.push_back( (size_t)segmentation_dims[i_dim] );
	}
	
	size_t Nz = data_size[3];
	size_t Ny = data_size[2];
	size_t Nx = data_size[1];

	for(int i=0; i<this->contrast_filled_volumes_.number(); ++i)
	{
		auto ptr_img = (CFImage*) contrast_filled_volumes_.sptr_image_wrap(i)->ptr_image();
		ptr_img->resize(Nx, Ny, Nz, 1);
		uint16_t const current_contast = ptr_img->getContrast();

		// #pragma omp parallel
		for( size_t nz=0; nz<Nz; nz++)
		for( size_t ny=0; ny<Ny; ny++)
		for( size_t nx=0; nx<Nx; nx++)
		{
			size_t linear_index_access = (nz*Ny + ny)*Nx + nx;
			std::vector<complex_float_t> curr_voxel = contrast_vector[linear_index_access];
			ptr_img->operator()(nx, ny, nz, 0) = curr_voxel[ current_contast ];						
		}
		
	}
	contrast_filled_volumes_.reorient(*(tlm_.get_sptr_geometry()));
}

void MRContrastGenerator::map_parameters()
{
	parameter_filled_volumes_.empty();
	parameter_filled_volumes_.push_back( get_parameter_map(0));
	parameter_filled_volumes_.push_back( get_parameter_map(1));
	parameter_filled_volumes_.push_back( get_parameter_map(2));
	parameter_filled_volumes_.push_back( get_parameter_map(3));
}


sirf::NiftiImageData3D<float> MRContrastGenerator::get_parameter_map(const int which_parameter) 
{
	this->tlm_.assign_tissues_to_labels();
	TissueVector tissue_params = this->tlm_.get_segmentation_tissues();

	LabelVolume parameter_map = this->tlm_.get_segmentation_labels();
	
	const int num_voxels = parameter_map.get_num_voxels();
	if(tissue_params.size() != num_voxels)
		throw std::runtime_error("Your Tissue label mapper gives a different number of voxels compared to your segmentation.");

	// #pragma omp parallel
	for( size_t i_vox=0; i_vox<num_voxels; i_vox++)
	{
		TissueParameter param_in_voxel = *(tissue_params[i_vox]);
		parameter_map(i_vox) = get_parameter_from_tissue(param_in_voxel, which_parameter);
	}

	return parameter_map;
}

float MRContrastGenerator::get_parameter_from_tissue(const TissueParameter tp, const int which_parameter) const
{
	switch(which_parameter)
	{
	case 0:
		return tp.mr_tissue_.spin_density_percentH2O_;
	case 1:
		return tp.mr_tissue_.t1_miliseconds_;
	case 2:
		return tp.mr_tissue_.t2_miliseconds_;
	case 3:
		return tp.mr_tissue_.cs_ppm_;
	default:
		throw std::runtime_error("Please as for parameter 0, 1, 2 or 3 and nothing else.");
	}
}

std::vector<complex_float_t> MRContrastGenerator::build_label_signal_map(std::vector<ExternalTissueSignal> ext_sig) const {

	TissueParameterList tpl = tlm_.get_tissue_parameter_list();

	if(ext_sig.size()<tpl.size())
		throw std::runtime_error("The external signal you supplied contains less labels than in the segmentation. Please supply one signal point for each label.");

	unsigned int max_label = 0;
	for(int i=0; i<tpl.size(); ++i)
		max_label = tpl[i].label_>max_label?tpl[i].label_:max_label;

	std::vector<complex_float_t> signal_map(max_label+1);

	for(int i=0; i<tpl.size(); ++i)
	{
		bool label_found = false;
		LabelType const curr_label = tpl[i].label_;
		for(int j=0; j<ext_sig.size(); ++j)
		{
			if(std::get<0>(ext_sig[j]) == curr_label)
			{
				signal_map[curr_label] = std::get<2>(ext_sig[j]);
				ext_sig.erase(ext_sig.begin() + j);
				label_found = true;
				break;
			}
		}
		if(!label_found)
			throw std::runtime_error("Could not find every label in the segmentation in the external signal dictionary.");
	}

	return signal_map;
}
void MRContrastGenerator::map_contrast(const std::vector<ExternalTissueSignal>& ext_sig)
{
	std::vector<complex_float_t> map_label_to_signal = build_label_signal_map(ext_sig);

	this->tlm_.assign_tissues_to_labels();
	this->contrast_filled_volumes_.empty();
	if(sptr_acqu_ == nullptr)
		throw std::runtime_error("Your contrast template acquisition is not set.");
		
	this->contrast_filled_volumes_ = GadgetronImagesVector(*sptr_acqu_);

	if( contrast_filled_volumes_.number() != 1)
		throw std::runtime_error("For external contrast signal mapping please supply template rawdata yiedling only one single image (i.e. only one contrast, repetition etc.)");

	TissueVector tissue_params = this->tlm_.get_segmentation_tissues();
	size_t const num_voxels = tissue_params.size();	

	std::vector< complex_float_t> magnetisation;
	magnetisation.resize(num_voxels);
	
	// #pragma omp parallel
	for (size_t i= 0; i<num_voxels; i++)
		magnetisation[i] = map_label_to_signal[tissue_params[i]->label_];
		
	const int* segmentation_dims = this->tlm_.get_segmentation_dimensions();

	size_t Nz = segmentation_dims[3];
	size_t Ny = segmentation_dims[2];
	size_t Nx = segmentation_dims[1];

	auto ptr_img = (CFImage*) contrast_filled_volumes_.sptr_image_wrap(0)->ptr_image();
	ptr_img->resize(Nx, Ny, Nz, 1);
	memcpy(ptr_img->begin(), &magnetisation[0], magnetisation.size() * sizeof(complex_float_t));

	contrast_filled_volumes_.reorient(*(tlm_.get_sptr_geometry()));
}

complex_float_t MRContrastGenerator::get_signal_for_tissuelabel( size_t const label )
{
	auto tissue_list = this->tlm_.get_tissue_parameter_list();

	TissueParameter tp;
	for(size_t i=0; i<tissue_list.size(); i++)
	{
		if( tissue_list[i].label_ == label)
		{
			tp = tissue_list[i];
			break;
		}

		if( i == tissue_list.size() -1)
		{
			std::stringstream error_msg_stream;
			error_msg_stream << "The tissue label: " << label << " is not part of the segmentation you gave to the simulation. Hence no SNR can be computed";
			throw std::runtime_error( error_msg_stream.str() );
		}
	}

	std::vector < complex_float_t >	(*contrast_map_function)(std::shared_ptr<TissueParameter> const ptr_to_tiss_par, const ISMRMRD::IsmrmrdHeader& ismrmrd_hdr);

	ISMRMRD::SequenceParameters const sequ_par = *(this->hdr_.sequenceParameters);
	std::string const sequ_name = *(sequ_par.sequence_type);

	if(sequ_name.compare("Flash") == 0)
	{
		contrast_map_function = &map_flash_contrast;
	}
	else if (sequ_name.compare("Bssfp") == 0)
	{
		contrast_map_function = &map_bssfp_contrast;
	}
	else
	{
		std::stringstream error_msg_stream;
		error_msg_stream << "The header you read in requires a contrast which has not been implemented yet. ";
		error_msg_stream << "The demanded sequence type is: " << sequ_name << ". ";
		error_msg_stream << "Please give another rawdata header or write the contrast map yourself and add an else if to the map_contrast method.";
		throw std::runtime_error( error_msg_stream.str() );
	}

	auto signal_vector = contrast_map_function( std::make_shared<TissueParameter>(tp), this->hdr_);

	return signal_vector[0];

};


std::vector < complex_float_t > map_flash_contrast(std::shared_ptr<TissueParameter> const ptr_to_tiss_par, const ISMRMRD::IsmrmrdHeader& ismrmrd_hdr)
{
	
	SequenceParameters const sequ_par = *(ismrmrd_hdr.sequenceParameters); 
	AcquisitionSystemInformation const asi = *(ismrmrd_hdr.acquisitionSystemInformation);
	
	SeqParamType const TE = *(sequ_par.TE);
	SeqParamType TR;
	
	try
	{	
		TR = *(sequ_par.echo_spacing);

	}
	catch(const std::runtime_error &e)
	{
		std::cout << "Caught exception in map_flash_contrast." << std::endl;
		std::cout << e.what() <<std::endl;
		std::cout << "Echo spacing was not set in header file, taking TR value instead." <<std::endl;
		TR = *(sequ_par.TR);
	}



	if ( TR.size() > 1 )
		throw std::runtime_error(" More than one echo spacing was given. Please give only one in Flash contrast.");
	
	SeqParamType const flip_angle_deg = *(sequ_par.flipAngle_deg);
	
	float const field_strength_t = *(asi.systemFieldStrength_T);
	

	if (flip_angle_deg.size() > 1)
		throw std::runtime_error(" More than one flip angle was given. Please give only one in Flash contrast.");

	size_t const num_echoes = TE.size();

	float const spin_dens = ptr_to_tiss_par->mr_tissue_.spin_density_percentH2O_;
	float const T1_ms = ptr_to_tiss_par->mr_tissue_.t1_miliseconds_;
	float const T2_ms = ptr_to_tiss_par->mr_tissue_.t2_miliseconds_;
	float const cs_ppm = ptr_to_tiss_par->mr_tissue_.cs_ppm_;
 
	std::vector< complex_float_t > contrast;
	contrast.resize( num_echoes );

	complex_float_t const imag_unit(0,1);
	float const gyro = 42.58 * 2*M_PI;

	// signal forumla
	for( int i_echo = 0; i_echo<num_echoes; i_echo++)
	{
		contrast[i_echo] = 	spin_dens * (float)sin( M_PI/180 * flip_angle_deg[0]) 
						 *(float)(1 - exp(-TR[0]/T1_ms)) / (float)( 1 - exp(-TR[0]/T1_ms)*cos(M_PI/180*flip_angle_deg[0]) )
						 *(float)exp( -TE[i_echo]/T2_ms) * exp(imag_unit * TE[i_echo] * gyro/1000.f * field_strength_t * cs_ppm);
	}

	
	return contrast;
}


std::vector <complex_float_t > map_bssfp_contrast( std::shared_ptr<TissueParameter> const ptr_to_tiss_par,
												   const ISMRMRD::IsmrmrdHeader& ismrmrd_hdr)
{
	// Signal model based on Haacke/Brown - Magnetic Resonance Imaging, Ch. 18.2.1, Eq. 18.57 f
	// Frequency response assumed to be of the form as described in  'Hargreaves et al. "Fat‐suppressed steady‐state free precession imaging using phase detection.", MRM(2003)' 

	SequenceParameters const sequ_par = *(ismrmrd_hdr.sequenceParameters);
	AcquisitionSystemInformation const asi = *(ismrmrd_hdr.acquisitionSystemInformation);

	SeqParamType const TE = *(sequ_par.TE);
	SeqParamType TR;

	try
	{	
		TR = *(sequ_par.echo_spacing);

	}
	catch(const std::runtime_error &e)
	{
		std::cout << "Caught exception in map_flash_contrast." << std::endl;
		std::cout << e.what() <<std::endl;
		std::cout << "Echo spacing was not set in header file, taking TR value instead." <<std::endl;
		TR = *(sequ_par.TR);
	}
	

	SeqParamType const flip_angle_deg = *(sequ_par.flipAngle_deg);
	float const field_strength_t = *(asi.systemFieldStrength_T);

	if (TR.size() > 1)
		throw std::runtime_error(" More than one echo spacing was given. Please give only one in Flash contrast.");

	if (flip_angle_deg.size() > 1)
		throw std::runtime_error(" More than one flip angle was given. Please give only one in Flash contrast.");

	size_t const num_echoes = TE.size();

	float const spin_dens = ptr_to_tiss_par->mr_tissue_.spin_density_percentH2O_;
	float const T1_ms = ptr_to_tiss_par->mr_tissue_.t1_miliseconds_;
	float const T2_ms = ptr_to_tiss_par->mr_tissue_.t2_miliseconds_;
	float const cs_ppm = ptr_to_tiss_par->mr_tissue_.cs_ppm_;

	std::vector< complex_float_t > contrast;
	contrast.resize( num_echoes );

	complex_float_t const imag_unit(0,1);
	float const gyro = 42.58*2*M_PI;

	float const E1 = exp( -1.f* TR[0]/T1_ms );
	float const E2 = exp( -1.f* TR[0]/T2_ms );

	float const sin_term = sin(M_PI/180.f*flip_angle_deg[0]);
	float const cos_term = cos(M_PI/180.f*flip_angle_deg[0]);
		

	// assuming passband step function behavior for magnetization sign.
	float const BSSFPPassbandWidth_Hz = 1000.f * 1.f/TR[0];	
	float const off_resonance_Hz = TR[0] * gyro/1000.f * field_strength_t * cs_ppm;
    bool const even_band = ((int)std::floor( (off_resonance_Hz + BSSFPPassbandWidth_Hz/2.f) / BSSFPPassbandWidth_Hz ) % 2) ==0;
    float const sign_bssfp = even_band? 1.f:-1.f;

	// assuming on resonance (or for fat, outside the signal-response valley)
	for( int i_echo = 0; i_echo<num_echoes; i_echo++)
	{
		contrast[i_echo] = sign_bssfp * spin_dens * sin_term * (1-E1) / ( 1-E1*cos_term - E2*(E1-cos_term) ) * (float)exp( -TE[i_echo]/T2_ms) * exp(imag_unit * TE[i_echo] * gyro/1000.f * field_strength_t * cs_ppm);
	}
	
	return contrast;
}



PETContrastGenerator::PETContrastGenerator (const LabelVolume& tissue_labels, const std::string& filename_tissue_parameter_xml) :
AbstractContrastGenerator(tissue_labels, filename_tissue_parameter_xml)
{
}

void PETContrastGenerator::set_template_image_from_file( const std::string& filename_header_with_ext )
{
 	 this->template_pet_image_data_ = sirf::STIRImageData(filename_header_with_ext);
 	 this->template_img_is_set_ = true;
 	 
};

std::vector< STIRImageData >& PETContrastGenerator::get_contrast_filled_volumes(bool const resample_output)
{
	if( resample_output == true )
		this->resample_to_template_image();

	return this->contrast_filled_volumes_;
}



std::vector< int > PETContrastGenerator::get_dimensions( void )
{
	if( this->template_img_is_set_ )
	{
		std::vector< int > dims;
		dims.resize(3, 0);

		auto works = this->template_pet_image_data_.get_dimensions(&dims[0]);

		if(works == -1)
			throw std::runtime_error("Irregular range of dimensions in PET image data.");

		std::reverse( dims.begin(), dims.end() );
		return dims;
	}
	else
	{	
		std::stringstream error_msg; error_msg << "From " << __FUNCTION__ << ": please set image template first using the dedicated function.";
		throw std::runtime_error( error_msg.str() );
	}
}

std::vector< float > PETContrastGenerator::get_voxel_sizes( void )
{
	if( this->template_img_is_set_ )
	{
		std::vector< float > vx_size;
		vx_size.resize(3, 0);

		this->template_pet_image_data_.get_voxel_sizes(&vx_size[0]);

		return vx_size;
	}
	else
	{	
		std::stringstream error_msg; error_msg << "From " << __FUNCTION__ << ": please set image template first using the dedicated function.";
		throw std::runtime_error( error_msg.str() );
	}
}

void PETContrastGenerator::map_tissue()
{
	this->map_contrast();
	this->map_attenuation();
}

void PETContrastGenerator::map_contrast()
{
	this->contrast_filled_volumes_.clear();
	this->map_tissueparams_member( CASE_MAP_PET_CONTRAST );
}

void PETContrastGenerator::map_attenuation()
{
	this->map_tissueparams_member( CASE_MAP_PET_ATTENUATION );
}

void PETContrastGenerator::map_tissueparams_member(int const case_map)
{
	using namespace stir;

	if (this->template_img_is_set_)
	{

		const int* segmentation_dims = this->tlm_.get_segmentation_dimensions();
		
		size_t const num_voxels = size_t( segmentation_dims[1] * segmentation_dims[2] * segmentation_dims[3] );

	  	std::vector < float > contrast_img;
	  	contrast_img.resize(num_voxels, 0);

	  	// STIRImageData pet_img_dat( template_pet_image_data_ );
	  	LabelVolume nifti_seg = this->tlm_.get_segmentation_labels();

		nifti_seg.write("temp_niftisegmentation.nii");
	  	STIRImageData pet_img_dat( "temp_niftisegmentation.nii" );

		std::vector< float > voxel_sizes(3,0.f);
		pet_img_dat.get_voxel_sizes(&voxel_sizes[0]);
		float const voxel_volume_ml = voxel_sizes[0] * voxel_sizes[1] * voxel_sizes[2] / 1000.f;

		TissueVector tissue_params = this->tlm_.get_segmentation_tissues();

		// #pragma omp parallel
		for( size_t i_vox=0; i_vox<num_voxels; i_vox++)
		{
			TissueParameter param_in_voxel = *(tissue_params[i_vox]);

			if(case_map==CASE_MAP_PET_CONTRAST)
				contrast_img[i_vox] = param_in_voxel.pet_tissue_.activity_kBq_ml_ * voxel_volume_ml ;						
			else if(case_map == CASE_MAP_PET_ATTENUATION)
				contrast_img[i_vox] = param_in_voxel.pet_tissue_.attenuation_1_by_cm_;										

		}
	
		pet_img_dat.set_data( &contrast_img[0] );
		this->contrast_filled_volumes_.push_back( pet_img_dat );
	}
	else
		throw std::runtime_error("To get dimensions of output correct please set image from file as template first using the dedicated method.");

}



void PETContrastGenerator::resample_to_template_image( void )
{

	NiftyResampler<float> resampler;

    resampler.set_interpolation_type_to_cubic_spline();
	resampler.set_reference_image(std::make_shared< sirf::STIRImageData > (this->template_pet_image_data_));
	
	auto contrast_volumes = this->get_contrast_filled_volumes();

	for(int i=0; i<contrast_volumes.size(); ++i)
	{
		resampler.set_floating_image(std::make_shared< sirf::STIRImageData >( contrast_volumes[i]));
		resampler.process();
		const std::shared_ptr<const sirf::ImageData> sptr_deformed_img = resampler.get_output_sptr();

		this->contrast_filled_volumes_[i] = sirf::STIRImageData( this->template_pet_image_data_ ); //constructor for STIRImageData from ImageData does not exist yet.
		
		sptr_deformed_img->copy(sptr_deformed_img->begin(),
								this->contrast_filled_volumes_[i].begin(), 
								this->contrast_filled_volumes_[i].end());
	}
}
