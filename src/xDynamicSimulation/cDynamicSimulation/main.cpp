/* ================================================

Author: Johannes Mayer
Date: 2018.03.15
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


// Main file to run test codes for DynamicsSimulation


#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>

#include "all_simulation_tests.h"

int recurse(int const idx, int const num_iter, std::vector< int > curr_perm, std::vector< size_t > dims)
{
	int l = (idx - curr_perm[num_iter])/dims[num_iter];
	return l;
}



int main( int argc, char *argv[] )
{

	try
	{
		if(argc > 1)
		{
			fprintf(stdout, "Please do not pass any arguments. This just runs test code.");
		}

		
		size_t const N=20;
		size_t const M=30;
		size_t const F=90;

		std::vector< size_t > dims;
		dims.push_back(N);
		dims.push_back(M);
		dims.push_back(F);

		std::vector< std::vector< int > > all_permuts;



		for( int i=0; i<N*M*F; i++)
		{
			std::vector< int > curr_permut;

			for( int j=0; j<dims.size(); j++)
			{
				int curr_idx = i;

				for(int k=0;k<j;k++)
				{
					// std::cout << "k " << k << " j "<< j << std::endl;
					curr_idx = recurse(curr_idx, k, curr_permut, dims);
				}

				curr_idx = curr_idx % dims[j];
				curr_permut.push_back(curr_idx);
			}
			all_permuts.push_back( curr_permut );
		}


		for(int i=0; i<all_permuts.size(); i++)
		{
			std::vector< int > curr_permut = all_permuts[i]; 
			for(int j=0; j<curr_permut.size(); j++)
				std::cout<< curr_permut[j]<< ",";
			std::cout<<std::endl;
		}


		// run_tests_auxiliary_testing_functions();
		// run_tests_auxiliary_input_output();
		// run_tests_tissueparameters();
		// run_tests_contrastgenerator();
		// run_tests_phantom_input();
		// run_tests_encoding();
		// run_tests_mr_acquisition_model();
		// run_tests_mr_dynamic_simulation();
		// run_tests_dynamics();

		return 0;
	}

	catch(const std::exception& e)
	{	
		std::cout << "Exception caught at highest level in main" << std::endl;
		std::cout<< e.what() << '\n';			
	}
	catch(...)
	{
		std::cout<< "An exception of unknown type was caught. The tests failed." <<std::endl;	
	}


}
