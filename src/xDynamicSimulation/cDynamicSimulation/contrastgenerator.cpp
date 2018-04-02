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
	this->rawdata_file_path = filepath_rawdata;
}

std::string AbstractContrastGenerator::get_rawdata_file_path( void )
{	
	if (this->rawdata_file_path.empty())
		throw std::runtime_error("Rawdata filepath is not set yet. Please do so first.");
	else
		return this -> rawdata_file_path;
}


MRContrastGenerator::MRContrastGenerator (LabelArray tissue_labels, std::string const filename_tissue_parameter_xml) :
AbstractContrastGenerator(tissue_labels, filename_tissue_parameter_xml)
{
	// constructor does already something
}


void MRContrastGenerator::read_rawdata_header()
{
	//Let's open the existing dataset
    ISMRMRD::Dataset d(this->rawdata_file_path.c_str(),"dataset", false);

    std::string xml;
    d.readHeader(xml);
    ISMRMRD::deserialize(xml.c_str(),this->hdr_);
	

	/*
    ISMRMRD::SequenceParameters sequ_par = this->hdr_.sequenceParameters.get();

	std::vector<float> TE = sequ_par.TE.get();

	size_t num_echoes = TE.size();
	
	for( int i=0; i<num_echoes; i++)
	{

		std::cout << "TE(" << i << ")" << TE[i] << std::endl;

	}*/    
}

void MRContrastGenerator::map_contrast()
{
	throw std::runtime_error(" this is not implemented yet, test it alter");
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
