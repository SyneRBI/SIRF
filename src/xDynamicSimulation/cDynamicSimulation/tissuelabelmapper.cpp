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



TissueArray assign_tissue_parameters_to_labels( TissueParameterList &tiss_list, LabelArray label_list )
{

	size_t num_tissue_params = tiss_list.size();

	std::map <int, TissueParameter* >  lut;

	for(int i =0; i<num_tissue_params; i++)
	{
		lut.insert(std::make_pair( tiss_list[i].label_, &tiss_list[i]));	//map label to pointer
	}

	TissueArray tiss_segm(boost::extents[1][1][1][1][1][1][1]);
	
	TissueParameter tiss_param;
	tiss_param.label_ = 0;


	return tiss_segm;

}


