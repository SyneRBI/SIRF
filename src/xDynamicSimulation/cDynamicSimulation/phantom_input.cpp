/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "sirf/cDynamicSimulation/phantom_input.h"

using namespace std;
using namespace H5;
using namespace sirf;


VoxelisedGeometricalInfo3D read_voxelised_geometry_info_from_h5_dataset( const std::string& h5_filename_with_suffix, const std::string& name_group )
{

	//  --------------------------- Extend to directions if necessary ---------------------------

	std::stringstream group_prefix;
	group_prefix << name_group << "/voxelised_geometry";

	std::stringstream namestream_size; 
	namestream_size << group_prefix.str() << "/size" ;

	std::stringstream namestream_offset; 
	namestream_offset << group_prefix.str() << "/offset" ;

	std::stringstream namestream_spacing; 
	namestream_spacing << group_prefix.str() << "/spacing"; 

	// ------------------------------------------------------------------------------------------

	H5T_class_t type_input_uint = H5T_INTEGER;
	PredType type_reader_uint = PredType::NATIVE_UINT32;

	std::vector< unsigned int > data_size = read_1D_dataset_from_h5< unsigned int >(h5_filename_with_suffix, namestream_size.str(), type_input_uint, type_reader_uint);

	H5T_class_t type_input_float = H5T_FLOAT;
	PredType type_reader_float = PredType::NATIVE_FLOAT;

	std::vector< float > data_offset  = read_1D_dataset_from_h5 <float> (h5_filename_with_suffix, namestream_offset.str(), type_input_float, type_reader_float );
	std::vector< float > data_spacing = read_1D_dataset_from_h5 <float> (h5_filename_with_suffix, namestream_spacing.str(), type_input_float, type_reader_float );


	if( data_size.size() !=3 || data_spacing.size() !=3 || data_offset.size() !=3)
		throw std::runtime_error( "The input is not three-dimensional geometry.");

	VoxelisedGeometricalInfo3D::Offset 	    geo_offset;
    VoxelisedGeometricalInfo3D::Spacing 	geo_spacing;
    VoxelisedGeometricalInfo3D::Size 		geo_size;
	
	for(int i=0;i<3;i++)
	{
		geo_offset [i] = data_offset[i];
		geo_spacing[i] = data_spacing[i];
		geo_size   [i] = data_size[i];
	}
	
	VoxelisedGeometricalInfo3D::Coordinate l_dir, p_dir, s_dir;

	l_dir[0]=1; l_dir[1]=0;	l_dir[2]=0;
	p_dir[0]=0; p_dir[1]=1;	p_dir[2]=0;
	s_dir[0]=0; s_dir[1]=0;	s_dir[2]=1;


    VoxelisedGeometricalInfo3D::DirectionMatrix geo_dir;
    geo_dir[0] = l_dir;    geo_dir[1] = p_dir;    geo_dir[2] = s_dir;

	VoxelisedGeometricalInfo3D geo_info(geo_offset, geo_spacing, geo_size, geo_dir);


	return geo_info;
}

sirf::NiftiImageData3D<float> read_nifti_from_h5( const std::string& h5_filename_with_suffix, const std::string& name_dataset, H5T_class_t data_type_dataset, H5::PredType data_type_reader )
{
	std::stringstream sstream_dataset;
	sstream_dataset << name_dataset <<  "/data";
	std::vector< float >	dat = read_1D_dataset_from_h5 <float> ( h5_filename_with_suffix, sstream_dataset.str(), data_type_dataset, data_type_reader );
	sirf::VoxelisedGeometricalInfo3D geo_info = read_voxelised_geometry_info_from_h5_dataset( h5_filename_with_suffix, name_dataset );

	sirf::NiftiImageData3D< float > nifti_img( &dat[0], geo_info);

	return nifti_img;
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

