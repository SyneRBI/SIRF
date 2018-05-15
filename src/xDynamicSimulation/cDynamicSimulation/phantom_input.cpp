/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "phantom_input.h"



using namespace H5;

ISMRMRD::NDArray< unsigned int > read_segmentation_from_h5( std::string const h5_filename_with_suffix)
{
	const H5std_string name_segmentation_dataset = "segmentation";

	H5File file( h5_filename_with_suffix, H5F_ACC_RDONLY );
	DataSet dataset = file.openDataSet( name_segmentation_dataset );

	H5T_class_t type_class = dataset.getTypeClass();

	if( type_class == H5T_INTEGER )
	{
		IntType intype = dataset.getIntType();

	    DataSpace dataspace = dataset.getSpace();

		hsize_t dimensions_input[ISMRMRD::ISMRMRD_NDARRAY_MAXDIM];
        hsize_t ndims = dataspace.getSimpleExtentDims( dimensions_input, NULL);

        for( int i = ndims; i < ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i++)
        	dimensions_input[i] = 1;

        std::vector < size_t > const input_dimensions (dimensions_input, dimensions_input + ISMRMRD::ISMRMRD_NDARRAY_MAXDIM);
	
		ISMRMRD::NDArray< unsigned int > segmentation( input_dimensions );

		dataset.read( segmentation.begin(), PredType::NATIVE_UINT32, dataspace, dataspace);

		return segmentation;
	}
	else
	{
		throw std::runtime_error("Please give an integer dataset as input for the segmentation.");
	}

}