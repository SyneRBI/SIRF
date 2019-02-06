/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "sirf/cDynamicSimulation/contrastgenerator.h"


#include <memory>
#include <stdexcept>
#include <cmath>
#include <math.h>
#include <algorithm>
//TODO Johannes #include <omp.h>

//#include "Testing/auxiliary_testing_functions.h"

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


MRContrastGenerator::MRContrastGenerator (const LabelVolume& tissue_labels, const std::string& filename_tissue_parameter_xml) :
AbstractContrastGenerator(tissue_labels, filename_tissue_parameter_xml)
{
}

void MRContrastGenerator::set_rawdata_header(const ISMRMRD::IsmrmrdHeader& hdr)
{
	this->hdr_ = hdr;
}

std::vector< ISMRMRD::Image< complex_float_t> >& MRContrastGenerator::get_contrast_filled_volumes()
{
	return this->contrast_filled_volumes_;	
}



void MRContrastGenerator::match_output_dims_to_headerinfo( void )
{
	

	using namespace ISMRMRD;

	std::vector< Image< complex_float_t > > padded_volumes;

    std::vector< Encoding > enc_vec = this->hdr_.encoding;

    if( enc_vec.size() != 1 )
    {
    	throw std::runtime_error(" More than one encoding object was given. ");
    }

    Encoding enc = enc_vec[0];
	EncodingSpace enc_space = enc.encodedSpace;
	MatrixSize enc_matrix_size = enc_space.matrixSize;	


	size_t num_contrast_volumes = this->contrast_filled_volumes_.size();



	for( int i_img=0; i_img<num_contrast_volumes; i_img++)
	{

		Image< complex_float_t > curr_image = this->contrast_filled_volumes_[i_img];
		auto padded_image = curr_image;
		padded_image.resize(enc_matrix_size.x, enc_matrix_size.y, enc_matrix_size.z, 1);
		


		for( size_t i_vox=0; i_vox< padded_image.getNumberOfDataElements();i_vox++)
			*(padded_image.begin() + i_vox) = std::complex< float > (0,0);

		std::vector < size_t > size_ratio; 
		size_ratio.push_back( enc_matrix_size.x / curr_image.getMatrixSizeX() );		
		size_ratio.push_back( enc_matrix_size.y / curr_image.getMatrixSizeY() );		
		size_ratio.push_back( enc_matrix_size.z / curr_image.getMatrixSizeZ() );		

		bool dims_match = true;
		for( int i = 0; i<3; i++)
			dims_match *= (size_ratio[i] == 1 || size_ratio[i] == 2);		
					
		if( dims_match )
		{
			for( size_t nz = 0; nz<curr_image.getMatrixSizeZ(); nz++ ) 
			for( size_t ny = 0; ny<curr_image.getMatrixSizeY(); ny++ ) 
			for( size_t nx = 0; nx<curr_image.getMatrixSizeX(); nx++ ) 
			{
				padded_image(nx + (size_ratio[0] - 1) * enc_matrix_size.x/4, ny + (size_ratio[1] - 1) * enc_matrix_size.y/4, nz + (size_ratio[3] - 1, 0) * enc_matrix_size.z/4, 0) = curr_image(nx, ny, nz, 0);
			}
		}
		else
		{
			throw std::runtime_error("The dimensions of the segmentation do not match the header information. Please modify either one to match the other. Dimension sizes in the segmentation half of the one in the encoded space are padded to cope for readout oversampling.");				
		}

		padded_volumes.push_back(padded_image);
	
	}
	this->contrast_filled_volumes_ = padded_volumes;
}


void MRContrastGenerator::map_contrast()
{
	this->tlm_.assign_tissues_to_labels();
	
	if( true ) 
	{
		std::vector< ISMRMRD::Image< complex_float_t> > temp(0);
		this->contrast_filled_volumes_.swap(temp);
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



	TissueVector tissue_params = this->tlm_.get_segmentation_tissues();
	size_t const num_voxels = tissue_params.size();	

	std::vector<std::vector< complex_float_t> > contrast_vector;
	contrast_vector.resize(num_voxels);

	
	#pragma omp parallel
	for (size_t i= 0; i<num_voxels; i++)
	{	
		contrast_vector[i] = contrast_map_function(tissue_params[i], this->hdr_);
		
	}

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


	ISMRMRD::Image< complex_float_t > contrast_img(Nx, Ny, Nz, 1);

	for( size_t i_contrast = 0; i_contrast<num_contrasts; i_contrast++)
	{
	
		// #pragma omp parallel
		for( size_t nz=0; nz<Nz; nz++)
		{
			for( size_t ny=0; ny<Ny; ny++)
			{
				for( size_t nx=0; nx<Nx; nx++)
				{
					size_t linear_index_access = (nz*Ny + ny)*Nx + nx;
					std::vector<complex_float_t> curr_voxel = contrast_vector[linear_index_access];
					contrast_img(nx, ny, nz, 0) = curr_voxel[i_contrast];						
				}
			}
		}

		contrast_img.setContrast(i_contrast);
		this->contrast_filled_volumes_.push_back( contrast_img );
	}
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
	float const gyro = 42.58;

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
	float const gyro = 42.58;

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

std::vector< STIRImageData >& PETContrastGenerator::get_contrast_filled_volumes()
{
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

		std::vector<size_t> data_dims;
		for( int i_dim=0; i_dim<3; i_dim++)
		{
			data_dims.push_back( (size_t)segmentation_dims[i_dim] );
		}
		
		TissueVector tissue_params = this->tlm_.get_segmentation_tissues();
	
		size_t Nz = data_dims[2];
		size_t Ny = data_dims[1];
		size_t Nx = data_dims[0];
		
		size_t const num_voxels = Nx*Ny*Nz;

	  	std::vector < float > contrast_img;
	  	contrast_img.resize(num_voxels, 0);


		STIRImageData pet_img_dat( template_pet_image_data_ );
		std::vector< float > voxel_sizes(3,0.f);
		pet_img_dat.get_voxel_sizes(&voxel_sizes[0]);
		float const voxel_volume_mm3 = voxel_sizes[0] * voxel_sizes[1] * voxel_sizes[2];


		// #pragma omp parallel
		for( size_t i_vox=0; i_vox<num_voxels; i_vox++)
		{
			TissueParameter param_in_voxel = *(tissue_params[i_vox]);

			if(case_map==CASE_MAP_PET_CONTRAST)
				contrast_img[i_vox] = (param_in_voxel.pet_tissue_.activity_kBq_ml_ * voxel_volume_mm3 / 1000.f);						
			else if(case_map == CASE_MAP_PET_ATTENUATION)
				{
					contrast_img[i_vox] = param_in_voxel.pet_tissue_.attenuation_1_by_cm_;						
				}
		}
	
		pet_img_dat.set_data( &contrast_img[0] );

		this->contrast_filled_volumes_.push_back( pet_img_dat );
	}
	else
		throw std::runtime_error("To get dimensions of output correct please set image from file as template first using the dedicated method.");

}


std::vector< float > PETContrastGenerator::get_template_based_volume_subset(const std::vector<float>& vol_data, const std::vector<size_t>& data_dims)
{
	std::vector< int > template_dims = this->get_dimensions();

	std::vector< float > out;
	out.resize(template_dims[0]*template_dims[1]*template_dims[2],0);

	std::vector< size_t > offsets;
	for(int i = 0; i<3; i++)
	{
		if(data_dims[i] >= template_dims[i])
			offsets.push_back(data_dims[i]/2 - template_dims[i]/2);
		else
			throw std::runtime_error("Please give only data which has equal or larger data dimensions than the template image.");
	}

	for(size_t nz = 0; nz<template_dims[2]; nz++)
	for(size_t ny = 0; ny<template_dims[1]; ny++)
	for(size_t nx = 0; nx<template_dims[0]; nx++)
	{
		
		size_t const linear_index_vol_data = ( (nz+offsets[2]) * data_dims[1] + (ny+offsets[1]) ) * data_dims[0] + (nx+offsets[0]);
		size_t const linear_index_subset = (nz*template_dims[1] + ny)*template_dims[0] + nx;
		
		out[linear_index_subset] = vol_data[linear_index_vol_data];
	}	

	return out;
}



