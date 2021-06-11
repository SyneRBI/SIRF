/* ================================================

Author: Johannes Mayer
Date: 2018.11.13
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once 


// for easier logging 
#define epiph(x) #x << " = " << x


#define SHARED_FOLDER_PATH "/media/sf_CCPPETMR/"
#define ANALYZE_OUTPUT_TESTPATH SHARED_FOLDER_PATH "analyze_test_output"

#define TESTDATA_PREFIX "TestData/Input/xDynamicSimulation/cDynamicSimulation/"
#define TESTDATA_OUT_PREFIX "TestData/Output/xDynamicSimulation/cDynamicSimulation/"

#define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH TESTDATA_PREFIX "Segmentations/xcat_phantom_incl_geomertry_128.h5"
#define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH TESTDATA_PREFIX "TemplateData/MR/CV_nav_cart_128Cube_FLASH_T1.h5"   
// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_nav_128_rpe_sfl_gc_usos8.h5" 

#define DISPLACEMENT_FIELD_PATH SHARED_FOLDER_PATH ""

#define PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "TemplateData/PET/template_image_input_contgen.hv"
#define PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "TemplateData/PET/template_image_input_acquisition.hv"
#define PET_TEMPLATE_ACQUISITION_DATA_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "TemplateData/PET/template_span11.hs"


#define TIME_POINTS_CARDIAC_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "SurrogateSignals/card_time"
#define CARDIAC_SIGNAL_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "SurrogateSignals/card_signal"

#define TIME_POINTS_RESP_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "SurrogateSignals/resp_time"
#define RESP_SIGNAL_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "SurrogateSignals/resp_signal"

#define XML_TEST_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "Segmentations/test_TissueParameters_XML.xml" 
#define XML_XCAT_PATH SHARED_FOLDER_PATH TESTDATA_PREFIX "Segmentations/XCAT_TissueParameters_XML.xml" 

// #define H5_PHANTOM_TEST_PATH  SHARED_FOLDER_PATH "h5_testfile_cube_size3.h5"
#define H5_PHANTOM_TEST_PATH SHARED_FOLDER_PATH "testdata_inputoutput/xcat_phantom_incl_geomertry_64.h5"

#define ACQU_FILE_NAME  SHARED_FOLDER_PATH "acquisitions_file_fwd_test.h5"

#define FILENAME_MR_RPE_SIM SHARED_FOLDER_PATH "testoutput_mr_rpe_simulation.h5"
#define FILENAME_MR_CONTRAST_DYNSIM  SHARED_FOLDER_PATH "testoutput_mr_dynamic_contrast_simulation.h5"
#define FILENAME_MR_MOTION_DYNSIM SHARED_FOLDER_PATH "testoutput_mr_dynamic_motion_simulation.h5"
#define FILENAME_MR_MOTION_CONTRAST_DYNSIM SHARED_FOLDER_PATH "testoutput_mr_dynamic_motion_contrast_simulation.h5"

#define FILENAME_STATICSIM_PET SHARED_FOLDER_PATH "testoutput_pet_static_simulation.hs"
#define FILENAME_MR_DEFORM_TEST SHARED_FOLDER_PATH "output_deforming_tests/deformed_img"



