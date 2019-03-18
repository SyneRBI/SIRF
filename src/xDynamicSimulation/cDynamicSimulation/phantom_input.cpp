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
	


	// TODO: read this from file aswell instead of hard coding
	VoxelisedGeometricalInfo3D::Coordinate l_dir, p_dir, s_dir;

	l_dir[0]=1; 	l_dir[1]= 0;	l_dir[2]=0;
	p_dir[0]=0; 	p_dir[1]= 1;	p_dir[2]=0;
	s_dir[0]=0; 	s_dir[1]= 0;	s_dir[2]=-1;


    VoxelisedGeometricalInfo3D::DirectionMatrix geo_dir;
    geo_dir[0] = l_dir;    geo_dir[1] = p_dir;    geo_dir[2] = s_dir;

	VoxelisedGeometricalInfo3D geo_info(geo_offset, geo_spacing, geo_size, geo_dir);


	return geo_info;
}

sirf::NiftiImageData3D<float> read_segmentation_to_nifti_from_h5(const std::string& h5_filename_with_suffix)
{
	std::cout << "Reading segmentation from file: " << h5_filename_with_suffix <<std::endl;

	std::string const dataset_name = "/segmentation";
	H5T_class_t type_input = H5T_INTEGER;
	H5::PredType type_reader = H5::PredType::NATIVE_UINT32;

	sirf::NiftiImageData3D<float> segmentation_nifti =  read_nifti_from_h5<DataTypeSegmentation>( h5_filename_with_suffix, dataset_name, type_input, type_reader );

	// sirf::NiftiImageData3D<float> read_nifti_from_h5( const std::string& h5_filename_with_suffix, const std::string& name_dataset, H5T_class_t data_type_dataset, H5::PredType data_type_reader )


	return segmentation_nifti;
}

std::vector< sirf::NiftiImageData3DDisplacement <float> > read_motionfields_to_nifti_from_h5(const std::string& h5_filename_with_suffix, const std::string& motionfield_type)
{

	std::cout << "Reading motionfield from file: " << h5_filename_with_suffix << std::endl;

	VoxelisedGeometricalInfo3D geo_info = read_voxelised_geometry_info_from_h5_dataset(h5_filename_with_suffix, "/segmentation");

	VoxelisedGeometricalInfo3D::Size const data_size = geo_info.get_size();
	size_t const num_voxels = data_size[0] * data_size[1] * data_size[2];

	VoxelisedGeometricalInfo3D::Spacing const spacing = geo_info.get_spacing();
	for( int i=0; i<3; ++i)
		std::cout <<  "#################SPACING in " << i << " = " <<  spacing[i] <<std::endl;


	std::stringstream sstream_name_dataset;
	sstream_name_dataset << "/motionfields/" << motionfield_type << "/data";

	H5T_class_t type_input_float = H5T_FLOAT;
	PredType type_reader_float = PredType::NATIVE_FLOAT;

	std::vector< float > dvf_data = read_1D_dataset_from_h5 <float> ( h5_filename_with_suffix, sstream_name_dataset.str(), type_input_float, type_reader_float );
	
	scale_vector_data_to_geometry( dvf_data, geo_info );

	size_t const num_dimensions = 3;
	size_t const num_phases = dvf_data.size() / (num_voxels * num_dimensions); 


	std::vector< sirf::NiftiImageData3DDisplacement <float> > out;

	for( int i=0; i<num_phases; i++ )
	{
		sirf::NiftiImageData3DDisplacement<float> dvf_i( &dvf_data[i * num_dimensions*num_voxels], geo_info);
		out.push_back(dvf_i);
	}

	return out;
}

void scale_vector_data_to_geometry( std::vector<float> &vect_data, const VoxelisedGeometricalInfo3D& geo_info)
{
	const VoxelisedGeometricalInfo3D::Spacing 	  input_spacing = geo_info.get_spacing();
    const VoxelisedGeometricalInfo3D::Size        input_size    = geo_info.get_size();

    size_t const Nz = input_size[2];
	size_t const Ny = input_size[1];
	size_t const Nx = input_size[0];

   	size_t const num_voxels = Nx*Ny*Nz;
	
	size_t const num_dimensions = 3;
	size_t const Nt = vect_data.size() / (num_voxels * num_dimensions); 
	
	if( vect_data.size() != num_voxels * Nt * num_dimensions)
		throw std::runtime_error("The vectorfield data does not match the geometrical information you passed.");

	std::array<float, 3> direction_sign{ 1, 1, 1};


	for( size_t nt=0; nt<Nt; nt++)
    for( size_t nz=0; nz<Nz; nz++)
    for( size_t ny=0; ny<Ny; ny++)
	for( size_t nx=0; nx<Nx; nx++)
	{
		size_t const voxel_access = (((nt*Nz + nz)*Ny + ny)*Nx + nx) *num_dimensions;
		for(size_t nv=0; nv<num_dimensions; nv++)
		{
			size_t const linear_access 		= voxel_access + nv;
			size_t const reverted_access 	= voxel_access + 2 - nv;

			// vect_data[linear_access] = direction_sign[nv] * input_spacing[nv] * vect_data[reverted_access];
			vect_data[linear_access] = direction_sign[nv] * input_spacing[nv] * vect_data[linear_access];
			// vect_data[linear_access] = direction_sign[nv] * 1 * vect_data[linear_access];
		}
	}
}



std::vector< sirf::NiftiImageData3DDisplacement <float> > read_cardiac_motionfields_to_nifti_from_h5( const std::string& h5_filename_with_suffix )
{
	return read_motionfields_to_nifti_from_h5(h5_filename_with_suffix, "cardiac");

}
std::vector< sirf::NiftiImageData3DDisplacement <float> > read_respiratory_motionfields_to_nifti_from_h5( const std::string& h5_filename_with_suffix )
{
	return read_motionfields_to_nifti_from_h5(h5_filename_with_suffix, "respiratory");
}



// OLD STUFF, DELETE AS SOON AS UP AND RUNNING


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

