/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once


#include <iostream>
#include <stdio.h>

#include <ismrmrd/ismrmrd.h>

#include <string>
#include <cstring>
#include <vector>
#include <typeinfo>


#include "H5Cpp.h"




typedef unsigned int DataTypeSegmentation;
typedef float DataTypeMotionFields;

ISMRMRD::NDArray< DataTypeSegmentation > read_segmentation_from_h5( std::string const h5_filename_with_suffix);
ISMRMRD::NDArray< DataTypeMotionFields > read_motionfield_from_h5( std::string const h5_filename_with_suffix, std::string const name_motion_field_dataset);

ISMRMRD::NDArray< DataTypeMotionFields > read_cardiac_motionfield_from_h5(std::string const h5_filename_with_suffix);
ISMRMRD::NDArray< DataTypeMotionFields > read_respiratory_motionfield_from_h5(std::string const h5_filename_with_suffix);




template< typename T > 
ISMRMRD::NDArray< T > read_dataset(std::string const h5_filename_with_suffix, std::string const name_dataset, H5T_class_t data_type_dataset, H5::PredType data_type_reader)
{
	using namespace H5;

	const H5std_string name_segmentation_dataset = name_dataset;

	H5File file( h5_filename_with_suffix, H5F_ACC_RDONLY );
	DataSet dataset = file.openDataSet( name_segmentation_dataset );

	H5T_class_t type_class = dataset.getTypeClass();

	if( type_class == data_type_dataset )
	{
		IntType intype = dataset.getIntType();

	    DataSpace dataspace = dataset.getSpace();

		hsize_t unsorted_dimensions_input[ISMRMRD::ISMRMRD_NDARRAY_MAXDIM];
        hsize_t ndims = dataspace.getSimpleExtentDims( unsorted_dimensions_input, NULL);

        for( int i = ndims; i < ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i++)
        	unsorted_dimensions_input[i] = 1;

      	hsize_t dimensions_input[ISMRMRD::ISMRMRD_NDARRAY_MAXDIM];

        for( int i = 0; i < ndims; i++)
			dimensions_input[i] = unsorted_dimensions_input[2-i];	        

		for( int i = ndims; i < ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i++)
        	dimensions_input[i] = 1;
        
        for( int i = 0; i < 7; i++)
			std::cout << dimensions_input[i] <<std::endl;

        std::vector < size_t > const input_dimensions (dimensions_input, dimensions_input + ISMRMRD::ISMRMRD_NDARRAY_MAXDIM);
	
		ISMRMRD::NDArray< T > output( input_dimensions );

		dataset.read( output.begin(), data_type_reader, dataspace, dataspace);

		return output;
	}
	else
	{
		throw std::runtime_error("Please give an integer dataset as input for the segmentation.");
	}
};