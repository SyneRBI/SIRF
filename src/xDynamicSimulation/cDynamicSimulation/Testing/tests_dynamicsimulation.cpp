/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


// Files collecting all tests for one specific module.


#include <stdio.h>
#include <iostream>
#include <stdexcept>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "dynamicsimulation_x.h"
#include "auxiliary_testing_functions.h"
#include "phantom_input.h"

#include "tests_dynamicsimulation.h"




bool tests_mr_dynsim::test_constructor( void ) 
{
	try
	{
	
	
	MRContrastGenerator mr_cont = aux_test::get_mock_mr_contrast_generator();

	MRDynamicSimulation mr_dyn_sim( MRContrastGenerator mr_cont_gen );

	return true;


	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}


}

void tests_mr_dynsim::test_extract_src_information( void )
{
	try
	{
	
	MRContrastGenerator mr_cont_gen = aux_test::get_mock_mr_contrast_generator();

	MRDynamicSimulation mr_dyn_sim( mr_cont_gen );

	mr_dyn_sim.set_filename_mr_rawdata( ISMRMRD_H5_TEST_PATH );

	mr_dyn_sim.extract_src_information();


	ISMRMRD::IsmrmrdHeader hdr = mr_dyn_sim.get_ismrmrd_header();

	std::stringstream xml;
	serialize(hdr, xml);

	std::cout << xml.str() << std::endl;	


	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}

}

bool tests_mr_dynsim::test_simulate_dynamics( void )
{
	try
	{
	
	ISMRMRD::NDArray< unsigned int > segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );
	MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

	MRDynamicSimulation mr_dyn_sim( mr_cont_gen );

	mr_dyn_sim.set_filename_mr_rawdata( ISMRMRD_H5_TEST_PATH );

	mr_dyn_sim.simulate_dynamics();

	return false;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}


