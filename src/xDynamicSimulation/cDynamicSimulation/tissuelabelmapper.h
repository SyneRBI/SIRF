/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once


#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
#include <ismrmrd/ismrmrd.h>
#include <map>

#include <utility>
#include <memory>

#include "tissueparameters.h"


typedef std::vector< TissueParameter* > TissueVector;
typedef ISMRMRD::NDArray<unsigned int> LabelArray;

class TissueLabelMapper{

public:
	TissueLabelMapper();

	std::string get_filepath_tissue_parameter_xml( void );
	void set_filepath_tissue_parameter_xml(std::string const filepath_tissue_parameter_xml);
	
	//void map_


private:

	bool flag_xml_path_is_set_ = false;

	std::string filepath_tissue_parameter_xml_;

	TissueParameterList tissue_parameter_list_;
	
	LabelArray segmentation_labels_;
	TissueVector segmentation_tissues_;
	
};


// public methods for the class
TissueVector assign_tissue_parameters_to_labels( TissueParameterList &tiss_list, LabelArray label_list );				
