/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "contrastgenerator.h"


#include <stdexcept>
#include <math.h>
#include <omp.h>

//#include "Testing/auxiliary_testing_functions.h"

AbstractContrastGenerator::AbstractContrastGenerator(LabelArray tissue_labels, std::string const filename_tissue_parameter_xml)
{
	this->tlm_ = TissueLabelMapper( tissue_labels, filename_tissue_parameter_xml );
	tlm_.map_labels_to_tissue_from_xml();

}


std::vector< ISMRMRD::Image< complex_float_t> > AbstractContrastGenerator::get_contrast_filled_volumes()
{
	return this->contrast_filled_volumes_;	
}


MRContrastGenerator::MRContrastGenerator (LabelArray tissue_labels, std::string const filename_tissue_parameter_xml) :
AbstractContrastGenerator(tissue_labels, filename_tissue_parameter_xml)
{
}

void MRContrastGenerator::set_rawdata_header(ISMRMRD::IsmrmrdHeader hdr)
{
	this->hdr_ = hdr;
}


void MRContrastGenerator::map_contrast()
{

	std::vector < complex_float_t >	(*contrast_map_function)(TissueParameter const * const ptr_to_tiss_par, ISMRMRD::IsmrmrdHeader * ptr_to_header);

	ISMRMRD::SequenceParameters sequ_par = this->hdr_.sequenceParameters.get(); 
	std::string const sequ_name = sequ_par.sequence_type.get();

	if(sequ_name.compare("Flash") == 0)
	{
		contrast_map_function = &map_flash_contrast;
	}
	else
	{
		throw std::runtime_error("The header you read in requires a contrast which has not been implemented yet. Please give another header or write the contrast map and add an else if to the map_contrast method.");
	}



	TissueVector tissue_params = this->tlm_.get_segmentation_tissues();
	size_t const num_voxels = tissue_params.size();	


	std::vector<std::vector< complex_float_t> > contrast_vector;
	contrast_vector.resize(num_voxels);

	
	//#pragma omp parallel
	for (size_t i= 0; i<num_voxels; i++)
	{	
		contrast_vector[i] = contrast_map_function(tissue_params[i], &(this->hdr_));
		
	}
	for(int i=0;i<8;i++)
	{
		std::vector< complex_float_t> first_voxel(contrast_vector[i]);
	}

	size_t const num_contrasts = contrast_vector[0].size();

	const size_t* segmentation_dims = this->tlm_.get_segmentation_dimensions();

	std::vector<size_t> data_size;
	
	for( int i_dim=0; i_dim<ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i_dim++)
	{
		data_size.push_back( segmentation_dims[i_dim] );
	}
	
		
	size_t Nz = data_size[2];
	size_t Ny = data_size[1];
	size_t Nx = data_size[0];


	ISMRMRD::Image< complex_float_t > contrast_img(Nx, Ny, Nz, 1);

	// sort data into NDArray
	
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


std::vector < complex_float_t > map_flash_contrast
( TissueParameter const * const ptr_to_tiss_par, ISMRMRD::IsmrmrdHeader * ptr_to_header)
{
	using namespace ISMRMRD;

	SequenceParameters sequ_par = ptr_to_header->sequenceParameters.get(); 
	AcquisitionSystemInformation asi = ptr_to_header->acquisitionSystemInformation.get();

	SeqParamType TE = sequ_par.TE.get();
	SeqParamType TR = sequ_par.TR.get();
	SeqParamType flip_angle_deg = sequ_par.flipAngle_deg.get();
	float const field_strength_t = asi.systemFieldStrength_T.get();

	if (TR.size() > 1)
		throw std::runtime_error(" More than one TR was given. Please give only one in Flash contrast.");

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
