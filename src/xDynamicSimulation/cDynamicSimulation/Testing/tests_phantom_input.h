/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>
#include <stdio.h>
#include <iostream>
#include <typeinfo>


#include <ismrmrd/ismrmrd.h>

#include "phantom_input.h"


// checking for 
template <class T>
bool check_array_content( ISMRMRD::NDArray <T> input_array, T central_value)
{
	try
	{

	unsigned int const NDims = input_array.getNDim();

	size_t const * dimensions = input_array.getDims();

	for( int idim = 0; idim < NDims; idim++)
	{
		std::cout << "dimensions" <<"[" << idim << "] = "<<  dimensions[idim] <<std::endl;

		if( dimensions[idim] % 2 == 0)
		{	
			throw std::runtime_error("Please only test arrays with odd number of elements in each dimension.");
		}
	}

	bool content_is_correct = true;
	
	for (int nk = 0; nk < dimensions[6]; nk++)	
	for (int nl = 0; nl < dimensions[5]; nl++)
	for (int nm = 0; nm < dimensions[4]; nm++)
	for (int nn = 0; nn < dimensions[3]; nn++)
	for (int nz = 0; nz < dimensions[2]; nz++)
	for (int ny = 0; ny < dimensions[1]; ny++)
	for (int nx = 0; nx < dimensions[0]; nx++)
			{	
				std::cout<<input_array(nx, ny, nz, nn, nm, nl, nk)<<std::endl;
				//content_is_correct *= (input_array(nx, ny, nz, nn, nm, nl, nk) == 0);
			}
	
	content_is_correct *= (input_array(dimensions[0]/2, dimensions[1]/2, dimensions[2]/2, 
										dimensions[3]/2, dimensions[5]/2, dimensions[5]/2, dimensions[6]/2) == central_value);


	return content_is_correct;
	
	}

	catch(std::runtime_error const &e)
	{
		std::cout << e.what() << std::endl;
	}
	
} 




bool test_read_dataset_from_h5(std::string h5_filename_with_suffix);
bool test_read_h5_segmentation_correct_dims( std::string h5_filename_with_suffix);
bool test_read_h5_segmentation_correct_content( std::string h5_filename_with_suffix);



