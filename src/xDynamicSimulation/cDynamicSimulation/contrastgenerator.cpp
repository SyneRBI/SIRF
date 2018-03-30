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

}
