/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include "tests_phantom_input.h"

bool test_read_h5_segmentation_correct_dims( std::string h5_filename_with_suffix)
{
	
	ISMRMRD::NDArray< unsigned int > segmentation = read_segmentation_from_h5( h5_filename_with_suffix );
	

	const size_t* dimensions = segmentation.getDims();
		
	size_t const input_seg_size = 3;
	bool dimensions_are_correct = ( dimensions[0] == input_seg_size) 
								* ( dimensions[1] == input_seg_size) 
								* ( dimensions[2] == input_seg_size);
	

	return dimensions_are_correct;
}


