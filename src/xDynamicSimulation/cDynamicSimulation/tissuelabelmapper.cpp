/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "sirf/cDynamicSimulation/tissuelabelmapper.h"



TissueLabelMapper::TissueLabelMapper() {}

TissueLabelMapper::TissueLabelMapper(const LabelVolume& segmentation_labels, const std::string& filepath_tissue_parameter_xml)
{
	this->filepath_tissue_parameter_xml_ = filepath_tissue_parameter_xml;
	this->segmentation_labels_ = segmentation_labels;

}


std::string TissueLabelMapper::get_filepath_tissue_parameter_xml()
{
	return this->filepath_tissue_parameter_xml_;
}


LabelVolume TissueLabelMapper::get_segmentation_labels( void )
{
	return this->segmentation_labels_;
}

TissueParameterList TissueLabelMapper::get_tissue_parameter_list( void )
{
	return this->tissue_parameter_list_;
}

void TissueLabelMapper::map_labels_to_tissue_from_xml( void )
{
	this->tissue_parameter_list_ = read_TissueParameters_from_xml(filepath_tissue_parameter_xml_);
	this-> assign_tissues_to_labels();
}


const int* TissueLabelMapper::get_segmentation_dimensions( void )
{
	return this->segmentation_labels_.get_dimensions();		
}


void TissueLabelMapper::assign_tissues_to_labels( void )
{
	this->segmentation_tissues_ = assign_tissue_parameters_to_labels( this->tissue_parameter_list_, this->segmentation_labels_);
}
	

TissueVector assign_tissue_parameters_to_labels( const TissueParameterList& tiss_list, const LabelVolume& label_volume )
{

	size_t const num_tissue_params = tiss_list.size();

	std::map <LabelType, std::shared_ptr<TissueParameter> >  lut;

	for(size_t i =0; i<num_tissue_params; i++)
	{
		lut.insert(std::make_pair( tiss_list[i].label_, std::make_shared<TissueParameter>(tiss_list[i]) ));	//map label to pointer
	}

	size_t const num_voxels = label_volume.get_num_voxels();
	
	TissueVector tissue_vect(num_voxels);
	
	for( size_t i_vox =0; i_vox<num_voxels; i_vox++)
	{
		auto key_value_pair = lut.find( (LabelType)label_volume( i_vox ) );
		if( key_value_pair != lut.end())
		{	
			tissue_vect[i_vox] = key_value_pair->second;
		}
		else
		{	
			std::stringstream msg;
			msg << "The label " <<  (LabelType)label_volume(i_vox) << " in your label volume does not appear in the label list.";
			throw std::runtime_error(msg.str());
		}

	}

	return tissue_vect;
}



void TissueLabelMapper::replace_petmr_tissue_parameters( const LabelType&  label, const TissueParameter& replacement_tiss)
{
	size_t const num_tissues = tissue_parameter_list_.size();

	bool label_found = false;

	for( size_t i_tiss=0; i_tiss<num_tissues; i_tiss++)
	{
		TissueParameter curr_tiss = this->tissue_parameter_list_[i_tiss];

		if(curr_tiss.label_ == label)
		{
			curr_tiss.mr_tissue_ = replacement_tiss.mr_tissue_;
			curr_tiss.pet_tissue_ = replacement_tiss.pet_tissue_;

			this->tissue_parameter_list_[i_tiss] = curr_tiss;			
			label_found = true;
			
			break;
		}
	}
	if( !label_found )
		throw std::runtime_error("The label you tried to replace did not exist in the segmentation.");
	// else
	// {       
	// 	this->segmentation_tissues_ = assign_tissue_parameters_to_labels( this->tissue_parameter_list_, segmentation_labels_);
	// }


}


