/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "phantom_input.h"



using namespace H5;


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

