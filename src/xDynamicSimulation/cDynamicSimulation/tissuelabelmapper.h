/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>
#include <vector>
#include <ismrmrd/ismrmrd.h>
#include <boost/multi_array.hpp>
#include <map>

#include <utility>


#include "tissueparameters.h"


typedef boost::multi_array< TissueParameter*, ISMRMRD::ISMRMRD_NDARRAY_MAXDIM> TissueArray;
typedef ISMRMRD::NDArray<unsigned int> LabelArray;

class TissueLabelMapper{

public:
	TissueLabelMapper();

	std::string get_filepath_tissue_parameter_xml( void );
	void set_filepath_tissue_parameter_xml(std::string const filepath_tissue_parameter_xml);
	


private:

	std::string filepath_tissue_parameter_xml_;

	TissueParameterList tissue_parameter_list_;
	
	LabelArray segmentation_labels_;
	TissueArray segmentation_tissues_;
	
};


// public methods for the class
TissueArray assign_tissue_parameters_to_labels( TissueParameterList &tiss_list, LabelArray label_list );				
