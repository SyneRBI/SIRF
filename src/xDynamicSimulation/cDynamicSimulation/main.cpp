/* ================================================

Author: Johannes Mayer
Date: 2018.03.15
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


// Main file to run test codes for DynamicsSimulation


#include <stdio.h>
#include <iostream>


#include "dynamicsimulation_x.h"

int main( int argc, char *argv[] )
{

	try
	{
		if(argc > 1)
		{
			fprintf(stdout, "Please do not pass any arguments. This just runs test code.");
		}

		return 0;
	}
	catch(...)
		{
		std::cout<< "An exception of unknown type was caught" <<std::endl;	
		}
}
