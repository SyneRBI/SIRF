/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */
// file containing auxiliary functions

#pragma once

#include <string>
#include <sstream>
#include <fstream>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "tissueparameters.h"
#include "tissuelabelmapper.h"



// for easier logging 
#define epiph(x) #x << " = " << x


// strings
#define XML_TEST_PATH "Testing/TestData/test_TissueParameters_XML.xml" 
#define ISMRMRD_H5_TEST_PATH "Testing/TestData/test_data_ismrmrd.h5"
#define H5_PHANTOM_TEST_PATH "Testing/TestData/h5_testfile_cube_size3.h5"
#define H5_XCAT_PHANTOM_PATH "Testing/TestData/xcat_tissue_segmentation_uint64.h5"

namespace aux_test
{

	TissueParameterList get_mock_tissue_param_list( void );
	LabelArray get_mock_label_array( void );

	TissueParameter get_mock_tissue_parameter( void );
	PETTissueParameter get_mock_PET_tissue_parameter( void );
	MRTissueParameter get_mock_MR_tissue_parameter( void );

	ISMRMRD::IsmrmrdHeader get_mock_ismrmrd_header( void );
	ISMRMRD::AcquisitionSystemInformation get_mock_acquisition_system_information( void );
	ISMRMRD::SequenceParameters get_mock_sequence_parameters( void );
	
	template <typename T>
	void write_ndarray_to_binary(std::string const output_name_without_ext, ISMRMRD::NDArray<T> data_array)
	{	
		std::cout<< "Writing file " <<output_name_without_ext << std::endl;
		std::stringstream name_stream;
		name_stream << output_name_without_ext << "_";

		const size_t* data_dimension = data_array.getDims();

		name_stream<< data_dimension[0];

		for(int i=1; i<7; i++)
		{	
			if( data_dimension[i] > 1)
			{
				name_stream << "x" << data_dimension[i];
			}
		}
		name_stream << ".raw";

		size_t num_elements = data_array.getNumberOfElements();
		std::vector <float> buffer;
		buffer.reserve(num_elements);
		
		T* arr_start = data_array.begin();

		for( size_t i=0; i<num_elements; i++)
		{
			buffer.push_back( std::abs(*(arr_start + i)) );
		}

		std::ofstream out;
		out.open( name_stream.str().c_str(), std::ios::out | std::ios::binary);

		out.write( reinterpret_cast<char*> (buffer.data()), buffer.size()*sizeof(float));

		out.close();

		std::cout<< "Finished writing file " << name_stream.str() << std::endl;


	};


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
