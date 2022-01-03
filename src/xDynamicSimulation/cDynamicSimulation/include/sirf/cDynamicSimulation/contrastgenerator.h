/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>
#include <memory>
#include <map>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "sirf/cDynamicSimulation/tissueparameters.h"
#include "sirf/cDynamicSimulation/tissuelabelmapper.h"

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/STIR/stir_data_containers.h"

#include "sirf/STIR/stir_types.h"

// base class for contrast generators. inherit for different modalities.
// Reading the header is the same for each modality (hopefully!!!).
// otherwise the header reading can be just inherited and implemented
// as a virutal method


#define CASE_MAP_PET_CONTRAST 0
#define CASE_MAP_PET_ATTENUATION 1


typedef std::vector<float> SeqParamType;

class AbstractContrastGenerator {

public:
	AbstractContrastGenerator(){};
	AbstractContrastGenerator(const LabelVolume& tissue_labels, const std::string& filename_tissue_parameter_xml);
	

	// pure virtual since formats are very diff for pet and mri and ct
	virtual void map_contrast()=0;
	void replace_petmr_tissue_parameters(LabelType label, TissueParameter tiss_param);	
	TissueParameter get_petmr_tissue_parameter( LabelType label );		
protected:
	virtual void resample_to_template_image( void )=0;
	TissueLabelMapper tlm_;

};

typedef std::tuple< LabelType, std::uint32_t, complex_float_t> ExternalTissueSignal; // label-#,  acquisition time-stamp, MR signal


class MRContrastGenerator : public AbstractContrastGenerator {

public:
	
	MRContrastGenerator (const LabelVolume& tissue_labels, const std::string& filename_tissue_parameter_xml);

	void set_template_rawdata(const sirf::MRAcquisitionData& ad);
	void set_rawdata_header(const ISMRMRD::IsmrmrdHeader& hdr);
	
	
	void map_contrast();
	void map_contrast(const std::vector<ExternalTissueSignal>& ext_sig);
	
	complex_float_t get_signal_for_tissuelabel( size_t const label );
	
	sirf::GadgetronImagesVector& get_contrast_filled_volumes(bool const resample_output=false);
	void set_contrast_filled_volumes(const sirf::GadgetronImagesVector& img)
	{
		contrast_filled_volumes_ = img;
	}

private:

	std::vector<complex_float_t> build_label_signal_map(std::vector<ExternalTissueSignal> ext_sig) const;

	void resample_to_template_image( void );
	sirf::GadgetronImagesVector contrast_filled_volumes_;
	ISMRMRD::IsmrmrdHeader hdr_;

	std::shared_ptr<sirf::MRAcquisitionData> sptr_acqu_;

};


std::vector < complex_float_t > map_flash_contrast( std::shared_ptr<TissueParameter> const ptr_to_tiss_par, 
													const ISMRMRD::IsmrmrdHeader& ismrmrd_hdr);


std::vector <complex_float_t > map_bssfp_contrast( std::shared_ptr<TissueParameter> const ptr_to_tiss_par,
													const ISMRMRD::IsmrmrdHeader& ismrmrd_hdr);


class PETContrastGenerator : public AbstractContrastGenerator {

public:

	PETContrastGenerator():AbstractContrastGenerator() {};

	PETContrastGenerator(const LabelVolume& tissue_labels, const std::string& filename_tissue_parameter_xml);


	void set_template_image_from_file( const std::string& filename_header_with_ext ); 

	std::vector< int > get_dimensions( void );
	std::vector< float > get_voxel_sizes( void );

	std::vector< sirf::STIRImageData >& get_contrast_filled_volumes(bool const resample_output=false);

	void map_tissue();
	void map_contrast();
	void map_attenuation();

private:

	void resample_to_template_image( void );

	bool template_img_is_set_ = false;

	std::vector < sirf::STIRImageData > contrast_filled_volumes_;
	void map_tissueparams_member(int const case_map);
	
	sirf::STIRImageData template_pet_image_data_;

};
