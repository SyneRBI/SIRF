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


void MRContrastGenerator::read_rawdata_header( std::string filename_ismrmrd_h5_file_with_ext ){}

void MRContrastGenerator::map_contrast(){}
