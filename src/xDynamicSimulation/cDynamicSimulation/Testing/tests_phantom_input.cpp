/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include "tests_phantom_input.h"


bool test_segmentation_dimensions_are_correct( std::string h5_filename_with_suffix)
{

	ISMRMRD::NDArray< unsigned int > segmenation = read_segmentation_from_h5( h5_filename_with_suffix );

	size_t dimensions[ISMRMRD_NDARRAY_MAXDIM] = segmentation.getDims();
	
	bool dimensions_are_correct = ( dimensions[0] == input_seg_size) 
								* ( dimensions[1] == input_seg_size) 
								* ( dimensions[2] == input_seg_size);
	

	return dimensions_are_correct;
}


