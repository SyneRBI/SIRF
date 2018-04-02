/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>
#include <stdexcept>
#include <math.h>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include <ismrmrd/dataset.h>

#include "tissueparameters.h"
#include "tissuelabelmapper.h"


// base class for contrast generators. inherit for different modalities.
// Reading the header is the same for each modality (hopefully!!!).
// otherwise the header reading can be just inherited and implemented
// as a virutal method

using ISMRMRD::NDArray;
using ISMRMRD::IsmrmrdHeader;


typedef std::vector<float> SeqParamType;


class AbstractContrastGenerator {

public:


	AbstractContrastGenerator(LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);
	

	// pure virtual since formats are very diff for pet and mri and ct
	virtual void read_rawdata_header( void ) = 0; 
	virtual void map_contrast() = 0;

	virtual std::string get_rawdata_file_path();
	virtual void set_rawdata_file_path(std::string filepath_rawdata);

protected:

	std::string rawdata_file_path;

	ISMRMRD::NDArray< complex_float_t > contrast_filled_volume_;
	TissueLabelMapper tlm_;

};



class MRContrastGenerator : public AbstractContrastGenerator {

public:
	
	MRContrastGenerator (LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);

	void read_rawdata_header();
	void map_contrast();


private:


	void get_sequence_type();
	std::string sequence_type_;
	IsmrmrdHeader hdr_;

};


std::vector < complex_float_t > map_flash_contrast( TissueParameter const * const ptr_to_tiss_par, 
													ISMRMRD::SequenceParameters * ptr_to_sequ_par);