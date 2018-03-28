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
	flag_xml_path_is_set_ = true;

}

std::string TissueLabelMapper::get_filepath_tissue_parameter_xml()
{
	if( flag_xml_path_is_set_ )
		return filepath_tissue_parameter_xml_;
	else
		throw std::runtime_error("Please first set the xml path using the setter");
}



TissueVector assign_tissue_parameters_to_labels( TissueParameterList &tiss_list, LabelArray label_volume )
{

	size_t num_tissue_params = tiss_list.size();

	std::map <int, TissueParameter* >  lut;

	for(int i =0; i<num_tissue_params; i++)
	{
		lut.insert(std::make_pair( tiss_list[i].label_, &tiss_list[i]));	//map label to pointer
	}

	size_t const num_voxels = label_volume.getNumberOfElements();

	TissueVector tissue_volume;
	tissue_volume.resize(num_voxels);

	for( int i_vox =0; i_vox<num_voxels; i_vox++)
	{
		auto key_value_pair = lut.find(label_volume(i_vox));
		if( key_value_pair != lut.end())
		{	
			tissue_volume[i_vox] = key_value_pair->second;
		}
		else
		{
			throw std::runtime_error("The label in your label volume does not appear in the label list.");
		}

	}

	return tissue_volume;

}


