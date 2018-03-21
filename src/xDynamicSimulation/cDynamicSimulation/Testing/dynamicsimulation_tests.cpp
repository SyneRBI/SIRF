/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


// Files collecting all tests for one specific module.


#include <stdio.h>
#include <iostream>
#include <stdexcept>


#include "dynamicsimulation_tests.h"


void run_tests_tissueparameters(void)
{
	bool tests_successful = false;




	// call every test here
	tests_successful = test_allocate_MRTissueParameter_successful();
	tests_successful = test_allocate_PETTissueParameter_successful();
	tests_successful = test_allocate_TissueParameter_successful();


	if ( !tests_successful )
	{
		throw std::runtime_error( "The tissueparameters tests failed.");
	}
	else
	{
		std::cout<< "The tissueparameters tests succeeded." <<std::endl;
	}
}


