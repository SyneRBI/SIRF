/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "sirf/cDynamicSimulation/phantom_input.h"
#include "auxiliary_testing_functions.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "tests_phantom_input.h"


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


