/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "auxiliary_testing_functions.h"

#include "sirf/common/GeometricalInfo.h"

#include "sirf/cDynamicSimulation/auxiliary_input_output.h"
#include "sirf/cDynamicSimulation/phantom_input.h"



#include "tests_phantom_input.h"


using namespace sirf;

bool test_read_h5_segmentation_correct_dims( std::string h5_filename_with_suffix)
{
	
	ISMRMRD::NDArray< DataTypeSegmentation > segmentation = read_segmentation_from_h5( h5_filename_with_suffix );
	

	const size_t* dimensions = segmentation.getDims();
	
	size_t const input_seg_size = 64;
	bool dimensions_are_correct = ( dimensions[0] == input_seg_size) 
								* ( dimensions[1] == input_seg_size) 
								* ( dimensions[2] == input_seg_size);
	

	return dimensions_are_correct;
}

bool test_read_h5_segmentation_correct_content( std::string h5_filename_with_suffix)
{
	
	ISMRMRD::NDArray< DataTypeSegmentation > segmentation = read_segmentation_from_h5( h5_filename_with_suffix );
	
	return check_array_content<DataTypeSegmentation>( segmentation);
		
}

void test_read_1D_dataset_from_h5( std::string h5_filename_with_suffix)
{

	std::string const name_dataset = "/segmentation/data";

	std::cout << "Reading dataset " << name_dataset <<std::endl;

	H5T_class_t type_input = H5T_INTEGER;
	H5::PredType type_reader = H5::PredType::NATIVE_UINT32;

 	std::vector< DataTypeSegmentation > input = read_1D_dataset_from_h5<DataTypeSegmentation>(h5_filename_with_suffix, name_dataset, type_input, type_reader);

	std::string output_name_xcat_seg =std::string( SHARED_FOLDER_PATH ) + "test_output_xcat_seg_from_1D_dataset" ;
	data_io::write_raw<DataTypeSegmentation> (output_name_xcat_seg, &input[0], input.size()); 	

}

void test_read_geometrical_info_from_h5( std::string h5_filename_with_suffix )
{
	std::string const group_name = "/segmentation";
 
	VoxelisedGeometricalInfo3D geom_info = read_voxelised_geometry_info_from_h5_dataset( h5_filename_with_suffix, group_name);

	VoxelisedGeometricalInfo3D::Offset 				input_offset 	 = geom_info.get_offset();
    VoxelisedGeometricalInfo3D::Spacing 			input_spacing	 = geom_info.get_spacing();
    VoxelisedGeometricalInfo3D::Size 				input_size   	 = geom_info.get_size() ;
    VoxelisedGeometricalInfo3D::DirectionMatrix 	input_direction  = geom_info.get_direction() ;


    for(int i=0; i<3; i++)
    {
    	std::cout << epiph(		input_offset[i]	    ) << std::endl;
    	std::cout << epiph(		input_spacing[i]	) << std::endl;
    	std::cout << epiph(		input_size[i]	 	) << std::endl;
    }

    std::cout << "Direction Matrix:" << std::endl;
    
    for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
    {
    	std::cout << epiph(		input_direction[i][j]	    ) << "   " ;
    	if (j == 2 )
    		std::cout << std::endl;
    }
}

void test_read_segmentation_to_nifti( std::string h5_filename_with_suffix )
{

	std::string const dataset_name = "/segmentation";
	H5T_class_t type_input = H5T_INTEGER;
	H5::PredType type_reader = H5::PredType::NATIVE_UINT32;

	// sirf::NiftiImageData3D<float> segmentation_nifti =  read_nifti_from_h5<DataTypeSegmentation>( h5_filename_with_suffix, dataset_name, type_input, type_reader );
	sirf::NiftiImageData3D<float> segmentation_nifti =  read_nifti_from_h5<DataTypeSegmentation>( h5_filename_with_suffix, dataset_name, type_input, type_reader );

	std::string output_name_seg_nifti =std::string( SHARED_FOLDER_PATH ) + "test_output_xcat_seg_from_nifti" ;
	
	std::cout << epiph( segmentation_nifti.get_max() ) <<std::endl;

	segmentation_nifti.write( output_name_seg_nifti);

}

void test_read_motionfield_to_nifti(  std::string h5_filename_with_suffix )
{

	std::string const type_motionfield = "cardiac";
	std::vector< sirf::NiftiImageData3DDisplacement <float> > all_dvfs = read_motionfields_to_nifti_from_h5(h5_filename_with_suffix, type_motionfield);

	std::cout << "Number of motion fields " << all_dvfs.size() << std::endl;

}























void test_read_h5_segmentation_for_xcat_input_check( std::string h5_filename_xcat_seg_with_suffix)
{
	ISMRMRD::NDArray< DataTypeSegmentation > segmentation = read_segmentation_from_h5(h5_filename_xcat_seg_with_suffix);

	std::string output_name_xcat_seg =std::string( SHARED_FOLDER_PATH ) + "test_output_xcat_seg_input_check" ;
	data_io::write_raw<DataTypeSegmentation> (output_name_xcat_seg, segmentation.begin(), segmentation.getNumberOfElements());

}


bool test_read_h5_motionfields( void )
{

try
	{

		ISMRMRD::NDArray< DataTypeMotionFields > resp_mvfs = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		ISMRMRD::NDArray< DataTypeMotionFields > card_mvfs = read_cardiac_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		
		auto resp_dims = resp_mvfs.getDims();

		for(int i=0; i<7; i++)
			std::cout << resp_dims[i] << std::endl;

		return true;
	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
	
}


