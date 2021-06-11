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


#include "sirf/common/GeometricalInfo.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"

#include "H5Cpp.h"




typedef unsigned int DataTypeSegmentation;
typedef float DataTypeMotionFields;


template < typename T >
std::vector< T > read_1D_dataset_from_h5( const std::string& h5_filename_with_suffix, const std::string& name_dataset, H5T_class_t data_type_dataset, H5::PredType data_type_reader  )
{
	using namespace H5;
	const H5std_string dataset_name_h5 = name_dataset;	

	std::cout << "Opening dataset " << name_dataset <<std::endl;

	H5File file( h5_filename_with_suffix, H5F_ACC_RDONLY );
	DataSet dataset = file.openDataSet( H5std_string( name_dataset ));

	H5T_class_t type_class = dataset.getTypeClass();

	if( type_class == data_type_dataset )
	{
	    DataSpace dataspace = dataset.getSpace();
		hsize_t dimensions_input[8];
        hsize_t ndims = dataspace.getSimpleExtentDims( dimensions_input, NULL);

        if( ndims != 1 )
        	throw std::runtime_error(" Please only read 1D data. Remaining information is stored in the geometry.");

        size_t const num_elements = dataspace.getSimpleExtentNpoints();

        std::cout << "Reading " << num_elements << " data elements." << std::endl;

        std::vector< T > output( num_elements );
		dataset.read( &output[0], data_type_reader, dataspace, dataspace);

		return output;
	}
	else
	{
		throw std::runtime_error("Please give read only from datasets with type passed to data_type.");
	}

}

sirf::VoxelisedGeometricalInfo3D read_voxelised_geometry_info_from_h5_dataset( const std::string& h5_filename_with_suffix, const std::string& name_group );

template< typename inputType >
sirf::NiftiImageData3D<float> read_nifti_from_h5( const std::string& h5_filename_with_suffix, const std::string& name_dataset, H5T_class_t data_type_dataset, H5::PredType data_type_reader )
{
	std::stringstream sstream_dataset;
	sstream_dataset << name_dataset <<  "/data";
	std::vector< inputType >	dat = read_1D_dataset_from_h5 <inputType> ( h5_filename_with_suffix, sstream_dataset.str(), data_type_dataset, data_type_reader );
	sirf::VoxelisedGeometricalInfo3D geo_info = read_voxelised_geometry_info_from_h5_dataset( h5_filename_with_suffix, name_dataset );

	sirf::NiftiImageData3D< float > nifti_img( &dat[0], geo_info);

	return nifti_img;
}

sirf::NiftiImageData3D<float> read_segmentation_to_nifti_from_h5(const std::string& h5_filename_with_suffix);
std::vector< sirf::NiftiImageData3DDisplacement <float> > read_motionfields_to_nifti_from_h5(const std::string& h5_filename_with_suffix, const std::string& motionfield_type);

void scale_vector_data_to_geometry( sirf::NiftiImageData3DDisplacement <float> &dvf );

std::vector< sirf::NiftiImageData3DDisplacement <float> > read_cardiac_motionfields_to_nifti_from_h5( const std::string& h5_filename_with_suffix );
std::vector< sirf::NiftiImageData3DDisplacement <float> > read_respiratory_motionfields_to_nifti_from_h5( const std::string& h5_filename_with_suffix );
