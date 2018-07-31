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


#include "stir_types.h"

// base class for contrast generators. inherit for different modalities.
// Reading the header is the same for each modality (hopefully!!!).
// otherwise the header reading can be just inherited and implemented
// as a virutal method

using ISMRMRD::NDArray;
using ISMRMRD::IsmrmrdHeader;


#define CASE_MAP_PET_CONTRAST 0
#define CASE_MAP_PET_ATTENUATION 1


typedef std::vector<float> SeqParamType;

class AbstractContrastGenerator {

public:


	AbstractContrastGenerator(LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);
	

	// pure virtual since formats are very diff for pet and mri and ct
	virtual void map_contrast()=0;

	

protected:

	TissueLabelMapper tlm_;

};



class MRContrastGenerator : public AbstractContrastGenerator {

public:
	
	MRContrastGenerator (LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);

	void set_rawdata_header(IsmrmrdHeader hdr);
	void map_contrast();

	virtual std::vector< ISMRMRD::Image< complex_float_t> > get_contrast_filled_volumes();

	void match_output_dims_to_headerinfo( void );

private:

	std::vector< ISMRMRD::Image< complex_float_t> > contrast_filled_volumes_;
	IsmrmrdHeader hdr_;

};


std::vector < complex_float_t > map_flash_contrast( TissueParameter const * const ptr_to_tiss_par, 
													ISMRMRD::IsmrmrdHeader * ptr_to_header);





#define CASE_MAP_CONTRAST 0
#define CASE_MAP_ATTENUATION 1


class PETContrastGenerator : public AbstractContrastGenerator {

public:

	PETContrastGenerator ( LabelArray tissue_labels, std::string const filename_tissue_parameter_xml );
	void set_rawdata_header ( void ) {};

	void map_contrast();
	void map_attenuation();

private:
	std::vector < Image3DF > contrast_filled_volumes_;
	map_tissueparams_member(int const case_map);
	// someheader;

}
