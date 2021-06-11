/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include <iostream>

#include "auxiliary_testing_functions.h"
#include "sirf/common/GeometricalInfo.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"
#include "sirf/cDynamicSimulation/phantom_input.h"
#include "tests_phantom_input.h"


using namespace sirf;
using namespace std;


bool test_read_1D_dataset_from_h5( std::string h5_filename_with_suffix)
{
	std::cout << "--- Running " << __FUNCTION__ <<std::endl;

	std::string const name_dataset = "/segmentation/data";

	std::cout << "Reading dataset " << name_dataset <<std::endl;

	H5T_class_t type_input = H5T_INTEGER;
	H5::PredType type_reader = H5::PredType::NATIVE_UINT32;

 	std::vector< DataTypeSegmentation > input = read_1D_dataset_from_h5<DataTypeSegmentation>(h5_filename_with_suffix, name_dataset, type_input, type_reader);

	stringstream ss_outname;
	ss_outname << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__;

	data_io::write_raw<DataTypeSegmentation> (ss_outname.str(), &input[0], input.size()); 	

	return true;
}

bool test_read_geometrical_info_from_h5( std::string h5_filename_with_suffix )
{
	std::cout << "--- Running " << __FUNCTION__ <<std::endl;

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
	return true;
}

bool test_read_segmentation_to_nifti( std::string h5_filename_with_suffix )
{
	std::cout << "--- Running " << __FUNCTION__ <<std::endl;

	std::string const dataset_name = "/segmentation";
	H5T_class_t type_input = H5T_INTEGER;
	H5::PredType type_reader = H5::PredType::NATIVE_UINT32;

	
	sirf::NiftiImageData3D<float> segmentation_nifti =  read_nifti_from_h5<DataTypeSegmentation>( h5_filename_with_suffix, dataset_name, type_input, type_reader );

	std::cout <<"Maximum val. in segmentation: " << epiph( segmentation_nifti.get_max() ) <<std::endl;

	stringstream ss_outname;
	ss_outname << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__;

	segmentation_nifti.write( ss_outname.str());

	return true;
}

bool test_read_motionfield_to_nifti(  std::string h5_filename_with_suffix )
{
	std::cout << "--- Running " << __FUNCTION__ <<std::endl;

	std::string const type_motionfield = "respiratory";
	std::vector< sirf::NiftiImageData3DDisplacement <float> > all_dvfs = read_motionfields_to_nifti_from_h5(h5_filename_with_suffix, type_motionfield);

	std::cout << "Number of motion fields " << all_dvfs.size() << std::endl;
	stringstream ss_outname;
	ss_outname << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__;


	all_dvfs[all_dvfs.size()-1].write( ss_outname.str() );

	return true;
}



