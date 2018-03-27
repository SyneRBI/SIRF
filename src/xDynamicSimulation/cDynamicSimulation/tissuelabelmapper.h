/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>
#include <vector>
#include <ismrmrd/ismrmrd.h>


#include "tissueparameters.h"


class TissueLabelMapper{

public:
	TissueLabelMapper();

	void set_filepath_tissue_parameter_xml(std::string const filepath_tissue_parameter_xml);
	void assign_tissue_parameters_to_labels( void );

private:

	std::string filepath_tissue_parameter_xml_;

	TissueParameterList tissue_parameter_list_;
	
	ISMRMRD::NDArray< unsigned int > segmentation_labels_;
	ISMRMRD::NDArray< TissueParameter* > segmentation_tissues_;
		
};
