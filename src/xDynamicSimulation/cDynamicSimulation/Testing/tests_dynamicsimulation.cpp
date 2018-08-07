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


using namespace sirf;



bool test_lin_combi_gen::test_get_all_combinations( void )
{

try
	{
		size_t const N = 1;
		size_t const M = 2;
		size_t const L = 3;

		std::vector< size_t > dims;
		dims.push_back(N);
		dims.push_back(M);
		dims.push_back(L);



		LinearCombiGenerator lcg( dims );

		lcg.compute_all_combinations();
		auto all_perm = lcg.get_all_combinations();

		for(size_t i=0;i<all_perm.size();i++)
		{
			auto curr_perm = all_perm[i];
			
			for(size_t j=0;j<curr_perm.size();j++)	
				std::cout << curr_perm[j] << "/";
			std::cout << std::endl;
		}


		bool test_succesful = (all_perm.size() == N*M*L);
		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}

}







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

	mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );

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

	mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );

	mr_dyn_sim.simulate_dynamics();

	mr_dyn_sim.write_simulation_results( FILENAME_DYNSIM );

	return true;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}


