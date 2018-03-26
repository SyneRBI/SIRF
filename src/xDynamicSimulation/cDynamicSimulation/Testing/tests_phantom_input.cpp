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

bool test_read_h5_segmentation_correct_content( std::string h5_filename_with_suffix)
{
	
	ISMRMRD::NDArray< unsigned int > segmentation = read_segmentation_from_h5( h5_filename_with_suffix );
	const size_t* dimensions = segmentation.getDims();

	int const center_position = dimensions[0] / 2 ;
	

	bool content_is_correct = true;

	for (int nz = 0; nz < dimensions[2]; nz++)
	for (int ny = 0; ny < dimensions[1]; ny++)
	for (int nx = 0; nx < dimensions[0]; nx++)
			{
				
				if( nx%3 == 0) { std::cout << "\n";}
				
				std::cout << segmentation(nx, ny, nz);

				if( nz==center_position && ny==center_position && nx==center_position)
				{
					content_is_correct *= (segmentation(nx, ny, nz) == 1);
				}
				else
				{
					content_is_correct *= (segmentation(nx, ny, nz) == 0);
				}
			}
	
	return content_is_correct;
}



