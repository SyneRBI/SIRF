/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tissuelabelmapper.h"



TissueLabelMapper::TissueLabelMapper() {}

TissueLabelMapper::TissueLabelMapper(LabelArray const segmentation_labels, std::string const filepath_tissue_parameter_xml)
{
	filepath_tissue_parameter_xml_ = filepath_tissue_parameter_xml;
	segmentation_labels_ = segmentation_labels;

}


std::string TissueLabelMapper::get_filepath_tissue_parameter_xml()
{
	return filepath_tissue_parameter_xml_;
}


LabelArray TissueLabelMapper::get_segmentation_labels( void )
{
	return segmentation_labels_;
}

void TissueLabelMapper::map_labels_to_tissue_from_xml( void )
{
	tissue_parameter_list_ = read_TissueParameters_from_xml(filepath_tissue_parameter_xml_);
	segmentation_tissues_ = assign_tissue_parameters_to_labels( tissue_parameter_list_, segmentation_labels_);
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


