/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "sirf/cDynamicSimulation/phantom_input.h"

#include <array>


using namespace H5;


sirf::VoxelisedGeometricalInfo3D read_voxelised_geometry_info_from_h5_dataset( const std::string& h5_filename_with_suffix, const std::string& name_group )
{

	//  --------------------------- Extend to directions if necessary ---------------------------

	std::stringstream namestream_size; 
	namestream_spacing << name_group << "/size" 

	std::stringstream namestream_offset; 
	namestream_spacing << name_group << "/offset" 

	std::stringstream namestream_spacing; 
	namestream_spacing << name_group << "/spacing" 

	// ---------------------------

	H5T_class_t type_input_uint = H5T_INTEGER;
	PredType type_reader_uint = PredType::NATIVE_UINT32;

	std::vector< unsigned int > data_size = read_1D_dataset_from_h5(h5_filename_with_suffix, namestream_size.str(), type_input_uint, type_reader_uint );


	H5T_class_t type_input_float = H5T_FLOAT;
	PredType type_reader_float = PredType::NATIVE_FLOAT;

	std::vector< float > data_offset  = read_1D_dataset_from_h5(h5_filename_with_suffix, namestream_offset.str(), type_input_float, type_reader_float );
	std::vector< float > data_spacing = read_1D_dataset_from_h5(h5_filename_with_suffix, namestream_spacing.str(), type_input_float, type_reader_float );

	for(int i =0; i<3; i++)
	{
		std::cout << data_size[i]  << std::endl;
		std::cout << data_offset[i] << std::endl;
		std::cout << data_spacing[i] << std::endl;
	}

	VoxelisedGeometricalInfo3D geo_info;

	return geo_info;
}




ISMRMRD::NDArray< DataTypeSegmentation > read_segmentation_from_h5( const std::string& h5_filename_with_suffix)
{
	std::string const name_dataset = "segmentation";

	std::cout << "Reading dataset /" << name_dataset <<std::endl;

	H5T_class_t type_input = H5T_INTEGER;
	PredType type_reader = PredType::NATIVE_UINT32;

 	return read_dataset< DataTypeSegmentation >(h5_filename_with_suffix, name_dataset, type_input, type_reader );
 
 }

ISMRMRD::NDArray< DataTypeMotionFields > read_motionfield_from_h5( const std::string& h5_filename_with_suffix, const std::string& name_motion_field_dataset )
{
	std::string const name_dataset =  "/motionfields/" + name_motion_field_dataset;

	std::cout << "Reading dataset /" << name_dataset <<std::endl;

	H5T_class_t type_input = H5T_FLOAT;
	PredType type_reader = PredType::NATIVE_FLOAT;

 	return read_dataset< DataTypeMotionFields >(h5_filename_with_suffix, name_dataset, type_input, type_reader );
}




ISMRMRD::NDArray< DataTypeMotionFields > read_cardiac_motionfield_from_h5( const std::string& h5_filename_with_suffix )
{
	std::string const name_motion_field_dataset = "cardiac";
	return read_motionfield_from_h5( h5_filename_with_suffix, name_motion_field_dataset);

}

ISMRMRD::NDArray< DataTypeMotionFields > read_respiratory_motionfield_from_h5( const std::string& h5_filename_with_suffix )
{
	std::string const name_motion_field_dataset = "respiratory";
	return read_motionfield_from_h5( h5_filename_with_suffix, name_motion_field_dataset);
}

