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


typedef std::vector<float> SeqParamType;


class AbstractContrastGenerator {

public:


	AbstractContrastGenerator(LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);
	

	// pure virtual since formats are very diff for pet and mri and ct
	virtual void map_contrast()=0;

	virtual ISMRMRD::NDArray< complex_float_t > get_contrast_filled_volume();

protected:

	ISMRMRD::NDArray< complex_float_t > contrast_filled_volume_;
	TissueLabelMapper tlm_;

};



class MRContrastGenerator : public AbstractContrastGenerator {

public:
	
	MRContrastGenerator (LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);

	void set_rawdata_header(ISMRMRD::IsmrmrdHeader hdr);
	void map_contrast();


private:

	IsmrmrdHeader hdr_;

};


std::vector < complex_float_t > map_flash_contrast( TissueParameter const * const ptr_to_tiss_par, 
													ISMRMRD::IsmrmrdHeader * ptr_to_header);
