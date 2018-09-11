/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */
// file containing auxiliary functions

#pragma once

#include <string>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "gadgetron_data_containers.h"
#include "gadgetron_image_wrap.h"

#include "auxiliary_input_output.h"

#include "phantom_input.h"
#include "tissueparameters.h"
#include "tissuelabelmapper.h"
#include "contrastgenerator.h"
#include "dynamics.h"
#include "encoding.h"

// for easier logging 
#define epiph(x) #x << " = " << x


// strings


// #define ISMRMRD_H5_TEST_PATH "/media/sf_SharedFolder/CCPPETMR/test_data_ismrmrd.h5"
// #define ISMRMRD_H5_TEST_PATH "/media/sf_SharedFolder/CCPPETMR/testdata_rpe128_ismrmrd.h5"

#define USE_64_CUBE_INPUT

#ifdef USE_64_CUBE_INPUT

	// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_SR_64Cube_1Echo_10Dyn.h5"
	#define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_nav_cart_64Cube_1Echo.h5"
	#define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_64_cubed.h5"
	#define DISPLACEMENT_FIELD_PATH SHARED_FOLDER_PATH "temp_folder_motion_dyn_0/motion_field_0.hdr"

#elif defined(USE_128_CUBE_INPUT)

	// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_nav_cart_128Cube_1Echo.h5"
	#define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_SR_128Cube_1Echo_3Dyn.h5"
	#define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_128_cubed.h5"
	#define DISPLACEMENT_FIELD_PATH SHARED_FOLDER_PATH ""

#endif



#define PET_TEMPLATE_IMAGE_DATA_PATH SHARED_FOLDER_PATH "pet_source_files/template_image_input.hv"
#define PET_TEMPLATE_ACQUISITION_DATA_PATH SHARED_FOLDER_PATH "pet_source_files/template_acquisition_input.hs"


#define XML_TEST_PATH SHARED_FOLDER_PATH "XMLTestData/test_TissueParameters_XML.xml" 
#define XML_XCAT_PATH SHARED_FOLDER_PATH "XMLTestData/XCAT_TissueParameters_XML.xml" 

// #define H5_PHANTOM_TEST_PATH  SHARED_FOLDER_PATH "h5_testfile_cube_size3.h5"
#define H5_PHANTOM_TEST_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_128_cubed.h5"


#define ACQU_FILE_NAME  SHARED_FOLDER_PATH "acquisitions_file_fwd_test.h5"

#define FILENAME_MR_RPE_SIM SHARED_FOLDER_PATH "testoutput_mr_rpe_simulation.h5"
#define FILENAME_MR_CONTRAST_DYNSIM  SHARED_FOLDER_PATH "testoutput_mr_dynamic_contrast_simulation.h5"
#define FILENAME_MR_MOTION_DYNSIM SHARED_FOLDER_PATH "testoutput_mr_dynamic_motion_simulation.h5"
#define FILENAME_MR_MOTION_CONTRAST_DYNSIM SHARED_FOLDER_PATH "testoutput_mr_dynamic_motion_contrast_simulation.h5"


#define FILENAME_DYNSIM_PET SHARED_FOLDER_PATH "testoutput_pet_dynamic_simulation.hs"


#define FILENAME_MR_DEFORM_TEST SHARED_FOLDER_PATH "output_deforming_tests/deformed_img"

// volume sizes

#define MOCK_FOV 256
#define MOCK_DATA_MATRIX_SIZE 128
#define MOCK_DATA_NUM_CHANNELS 1
#define MOCK_DATA_RO_OVERSAMPLING 2
#define MOCK_IMAGE_TYPE 5 // from ismrmrd enum ISMRMRD_IMTYPE_COMPLEX   = 5
#define MOCK_DATA_TYPE 7 // from ismrmrd enum ISMRMRD_CXFLOAT = 7
#define MOCK_FIELD_STRENGTH 1

// mock signal
#define MOCK_NUM_SIG_PTS 10


namespace aux_test
{

	TissueParameterList get_mock_tissue_param_list( void );
	LabelArray get_mock_label_array( void );

	TissueParameter get_mock_tissue_parameter( void );
	PETTissueParameter get_mock_PET_tissue_parameter( void );
	MRTissueParameter get_mock_MR_tissue_parameter( void );

	std::pair< TissueParameter, TissueParameter> get_mock_contrast_signal_extremes( void );

	MRContrastGenerator get_mock_mr_contrast_generator( void );
	PETContrastGenerator get_mock_pet_contrast_generator( void );

	ISMRMRD::IsmrmrdHeader get_mock_ismrmrd_header( void );
	std::string get_serialized_mock_ismrmrd_header( void );


	ISMRMRD::AcquisitionSystemInformation get_mock_acquisition_system_information( void );
	ISMRMRD::SequenceParameters get_mock_sequence_parameters( void );
	ISMRMRD::ExperimentalConditions get_mock_experimental_conditions( void );
	std::vector< ISMRMRD::Encoding > get_mock_encoding_vector( void );

	ISMRMRD::ExperimentalConditions get_mock_experimental_conditions( void );
	ISMRMRD::EncodingSpace get_mock_encoded_space( void );
	ISMRMRD::EncodingSpace get_mock_recon_space( void );
	ISMRMRD::EncodingLimits get_mock_encoding_limits( void );



	ISMRMRD::NDArray<complex_float_t> get_mock_ndarray_with_cube( void );
	ISMRMRD::Image< complex_float_t > get_mock_ismrmrd_image_with_cube( void );

	ISMRMRD::NDArray<complex_float_t> get_mock_csm( void );
	sirf::CoilDataAsCFImage get_mock_coildata_as_cfimage( void );

	ISMRMRD::AcquisitionHeader get_mock_acquisition_header( void );	
	sirf::AcquisitionsVector get_mock_acquisition_vector ( ISMRMRD::IsmrmrdHeader );	

	sirf::RPETrajectoryContainer get_mock_radial_trajectory(size_t const NRad, size_t const NAng);


	SignalContainer get_mock_motion_signal( void );
	SignalContainer get_mock_sinus_signal( sirf::AcquisitionsVector acq_vec);
	SignalContainer get_mock_contrast_signal( sirf::AcquisitionsVector acq_vec);

	



	template <typename T> bool equal_array_content( ISMRMRD::NDArray<T> one_array, ISMRMRD::NDArray<T> other_array)
	{
		size_t const num_elements = one_array.getNumberOfElements();
		size_t const num_elements_other = other_array.getNumberOfElements();

		if( num_elements != num_elements_other)
		{
			return false;
		}
		else
		{
			bool content_is_equal = true;
			for(int i=0; i<num_elements; i++)
			{
				content_is_equal *= ( one_array(i) == other_array(i) );
			}
			return content_is_equal;
		}
	};


}// namespace aux_test
