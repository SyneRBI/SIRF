/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "tissueparameters.h"
#include "tissuelabelmapper.h"


// base class for contrast generators. inherit for different modalities.
// Reading the header is the same for each modality (hopefully!!!).
// otherwise the header reading can be just inherited and implemented
// as a virutal method

using ISMRMRD::NDArray;
using ISMRMRD::IsmrmrdHeader;


class AbstractContrastGenerator {

public:


	AbstractContrastGenerator(LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);
	

	// pure virtual since formats are very diff for pet and mri and ct
	virtual void read_rawdata_header(std::string filename_ismrmrd_h5_file_with_ext ) = 0; 
	virtual void map_contrast() = 0;

	virtual std::string get_rawdata_file_path();
	virtual void set_rawdata_file_path(std::string filepath_rawdata);

protected:

	ISMRMRD::NDArray<float> contrast_filled_volume_;
	TissueLabelMapper tlm_;

};



class MRContrastGenerator : public AbstractContrastGenerator {

public:
	
	MRContrastGenerator (LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);

	void read_rawdata_header( std::string filename_ismrmrd_h5_file_with_ext );
	void map_contrast();


private:

	IsmrmrdHeader hdr_;

};