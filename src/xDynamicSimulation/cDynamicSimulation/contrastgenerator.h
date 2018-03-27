/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#pragma once



class ContrastGenerator{

public:
	ContrastGenerator();

	void set_filepath_tissue_parameter_xml_(std::string filepath_tissue_parameter_xml);

private:

	std::string filepath_tissue_parameter_xml_;

	std::vector < TissueParameter > tissue_parameter_list_;
	ISMRMRD::NDArray< unsigned int > segmentation_labels_;
	ISMRMRD::NDArray< TissueParameter* > ;


};
