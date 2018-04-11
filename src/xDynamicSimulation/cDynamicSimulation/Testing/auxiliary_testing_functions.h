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

#include "tissueparameters.h"
#include "tissuelabelmapper.h"



// for easier logging 
#define epiph(x) #x << " = " << x


// strings
#define XML_TEST_PATH "Testing/TestData/test_TissueParameters_XML.xml" 
#define XML_XCAT_PATH "Testing/TestData/XCAT_TissueParameters_XML.xml" 
#define ISMRMRD_H5_TEST_PATH "/media/sf_SharedFiles/test_data_ismrmrd.h5"
#define H5_PHANTOM_TEST_PATH "Testing/TestData/h5_testfile_cube_size3.h5"
#define H5_XCAT_PHANTOM_PATH "Testing/TestData/xcat_tissue_segmentation_uint64.h5"


// volume sizes

#define MOCK_DATA_MATRIX_SIZE 64
#define MOCK_DATA_NUM_CHANNELS 2
#define MOCK_DATA_RO_OVERSAMPLING 1
#define MOCK_IMAGE_TYPE 5 // from ismrmrd enum ISMRMRD_IMTYPE_COMPLEX   = 5



namespace aux_test
{

	TissueParameterList get_mock_tissue_param_list( void );
	LabelArray get_mock_label_array( void );

	TissueParameter get_mock_tissue_parameter( void );
	PETTissueParameter get_mock_PET_tissue_parameter( void );
	MRTissueParameter get_mock_MR_tissue_parameter( void );


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
	CoilDataAsCFImage get_mock_coildata_as_cfimage( void );

	ISMRMRD::AcquisitionHeader get_mock_acquisition_header( void );	
	ISMRMRD::Acquisition get_mock_ismrmrd_acquisition ( void );	


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
