/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "phantom_input.h"
#include "auxiliary_testing_functions.h"
#include "../auxiliary_input_output.h"

#include "tests_phantom_input.h"


bool test_read_h5_segmentation_correct_dims( std::string h5_filename_with_suffix)
{
	
	ISMRMRD::NDArray< unsigned int > segmentation = read_segmentation_from_h5( h5_filename_with_suffix );
	

	const size_t* dimensions = segmentation.getDims();
	
	size_t const input_seg_size = 33;
	bool dimensions_are_correct = ( dimensions[0] == input_seg_size) 
								* ( dimensions[1] == input_seg_size) 
								* ( dimensions[2] == input_seg_size);
	

	return dimensions_are_correct;
}

bool test_read_h5_segmentation_correct_content( std::string h5_filename_with_suffix)
{
	
	ISMRMRD::NDArray< unsigned int > segmentation = read_segmentation_from_h5( h5_filename_with_suffix );
	
	return check_array_content<unsigned int>( segmentation);
		
}


void test_read_h5_segmentation_for_xcat_input_check( std::string h5_filename_xcat_seg_with_suffix)
{
	ISMRMRD::NDArray< unsigned int > segmentation = read_segmentation_from_h5(h5_filename_xcat_seg_with_suffix);

	std::string output_name_xcat_seg = "/media/sf_SharedFiles/test_output_xcat_seg_input_check";
	data_io::write_raw<unsigned int> (output_name_xcat_seg, segmentation.begin(), segmentation.getNumberOfElements());

}


