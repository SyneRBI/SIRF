/* ================================================

Author: Johannes Mayer
Date: 2018.03.15
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


// Main file to run test codes for DynamicsSimulation


#include <stdio.h>
#include <iostream>


#include "all_simulation_tests.h"

int main( int argc, char *argv[] )
{

	try
	{
		if(argc > 1)
		{
			fprintf(stdout, "Please do not pass any arguments. This just runs test code.");
		}

		
		// run_tests_auxiliary_testing_functions();
		// run_tests_auxiliary_input_output();
		// run_tests_tissueparameters();
		// run_tests_contrastgenerator();
		run_tests_phantom_input();
		// run_tests_encoding();
		// run_tests_mr_acquisition_model();
		// run_tests_mr_dynamic_simulation();

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
