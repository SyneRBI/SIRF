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

#include "auxiliary_testing_functions.h"


// checking for 
template <class T>
bool check_array_content( ISMRMRD::NDArray <T> input_array)
{
	try
	{

	unsigned int const NDims = input_array.getNDim();

	size_t const * dimensions = input_array.getDims();

	for( int idim = 0; idim < NDims; idim++)
	{
		std::cout << epiph(dimensions[idim]) <<std::endl;

		if( dimensions[idim] % 2 == 0)
		{	
			throw std::runtime_error("Please only test arrays with odd number of elements in each dimension.");
		}
	}

	bool content_is_correct = true;
	

	size_t Nk = dimensions[6];
	size_t Nl = dimensions[5];
	size_t Nm = dimensions[4];
	size_t Nn = dimensions[3];
	size_t Nz = dimensions[2];
	size_t Ny = dimensions[1];
	size_t Nx = dimensions[0];





	for (int nk = 0; nk < Nk; nk++)	
	for (int nl = 0; nl < Nl; nl++)
	for (int nm = 0; nm < Nm; nm++)
	for (int nn = 0; nn < Nn; nn++)
	for (int nz = 0; nz < Nz; nz++)
	for (int ny = 0; ny < Ny; ny++)
	for (int nx = 0; nx < Nx; nx++)
			{	
				size_t current_access = (((((nk*Nl + nl)*Nm + nm)*Nn + nn)*Nz + nz)*Ny + ny)*Nx+nx;
				std::cout << epiph(current_access) << "   " << epiph(input_array(nx, ny, nz, nn, nm, nl, nk))<<std::endl;
				content_is_correct *= (input_array(nx,ny,nz,nn,nm,nl,nk) == current_access + 1);
				
			}
	

	return content_is_correct;
	
	}

	catch(std::runtime_error const &e)
	{
		std::cout << e.what() << std::endl;
	}
	
} 



void test_read_h5_segmentation_for_xcat_input_check(std::string h5_filename_with_suffix);
bool test_read_dataset_from_h5(std::string h5_filename_with_suffix);
bool test_read_h5_segmentation_correct_dims( std::string h5_filename_with_suffix);
bool test_read_h5_segmentation_correct_content( std::string h5_filename_with_suffix);



