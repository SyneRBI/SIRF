/*
author: johannes mayer
date: 15. March 2018

*/


#include "dynamicsimulation_x.h"



void foo( void )
{
	
	ImagesReconstructor ir;		
	const char* nameOfClass = ir.class_name();

	fprintf(stdout, "This function is called from the cdynamicsimulation library.\n");
	fprintf(stdout, "It created a class called %s from gadgetron_x.h.\n", nameOfClass);

}
