/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "contrastgenerator.h"


AbstractContrastGenerator::AbstractContrastGenerator(LabelArray tissue_labels, std::string const filename_tissue_parameter_xml)
{
	this->tlm_ = TissueLabelMapper( tissue_labels, filename_tissue_parameter_xml );
	tlm_.map_labels_to_tissue_from_xml();

}

void AbstractContrastGenerator::set_rawdata_file_path(std::string const filepath_rawdata)
{
	this->rawdata_file_path_ = filepath_rawdata;
}

std::string AbstractContrastGenerator::get_rawdata_file_path( void )
{	
	if (this->rawdata_file_path_.empty())
		throw std::runtime_error("Rawdata filepath is not set yet. Please do so first.");
	else
		return this -> rawdata_file_path_;
}


ISMRMRD::NDArray< complex_float_t > AbstractContrastGenerator::get_contrast_filled_volume()
{
	return this->contrast_filled_volume_;	
}


MRContrastGenerator::MRContrastGenerator (LabelArray tissue_labels, std::string const filename_tissue_parameter_xml) :
AbstractContrastGenerator(tissue_labels, filename_tissue_parameter_xml)
{
}


void MRContrastGenerator::read_rawdata_header()
{
	//Let's open the existing dataset
    ISMRMRD::Dataset d(this->rawdata_file_path_.c_str(),"dataset", false);

    std::string xml;
    d.readHeader(xml);
    ISMRMRD::deserialize(xml.c_str(),this->hdr_);

}

void MRContrastGenerator::map_contrast()
{

	read_rawdata_header();
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
	

	for (size_t i= 0; i<num_voxels; i++)
	{
		contrast_vector[i] = contrast_map_function(tissue_params[i], &(this->hdr_));
	}
	size_t const num_echoes = contrast_vector[0].size();

	const size_t* segmentation_dims = this->tlm_.get_segmentation_dimensions();

	std::vector<size_t> data_size;
	data_size.resize(ISMRMRD::ISMRMRD_NDARRAY_MAXDIM);
	for( int i_dim=0; i_dim<ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i_dim++)
	{
		data_size[i_dim] = segmentation_dims[i_dim];
	}

	data_size[3] = num_echoes;

	this->contrast_filled_volume_.resize(data_size);
	
	for( size_t i_echo = 0; i_echo<num_echoes; i_echo++)
		for

	

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
						   	*(float)exp( -TE[i_echo]/T2_ms) * exp(imag_unit * TE[i_echo] * gyro/1000.f * field_strength_t);
	}

	return contrast;
}
