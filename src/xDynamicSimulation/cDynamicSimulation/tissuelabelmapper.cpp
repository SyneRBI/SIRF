/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tissuelabelmapper.h"



TissueLabelMapper::TissueLabelMapper() {}

void TissueLabelMapper::set_filepath_tissue_parameter_xml(std::string const filepath_tissue_parameter_xml)
{
	filepath_tissue_parameter_xml_ = filepath_tissue_parameter_xml;
}

std::string TissueLabelMapper::get_filepath_tissue_parameter_xml()
{
	return filepath_tissue_parameter_xml_;
}



void TissueLabelMapper::assign_tissue_parameters_to_labels( void )
{/*
	tissue_parameter_list_ = read_TissueParameters_from_xml(filepath_tissue_parameter_xml_);
			
	typedef std::map <int, TissueParameter* > LabelTissueMap;

	size_t num_tissue_params = tissue_parameter_list_.size();

	LabelTissueMap lut;

	for(int i =0; i<num_tissue_params; i++)
	{
		lut.insert(std::make_pair( tissue_parameter_list_[i].label_, &tissue_parameter_list_[i]);	)	
	}
*/
}


