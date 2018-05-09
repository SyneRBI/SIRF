/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


// Files collecting all tests for one specific module.


#include <stdio.h>
#include <iostream>
#include <stdexcept>

#include "tests_auxiliary_testing_functions.h"
#include "tests_auxiliary_input_output.h"
#include "tests_tissueparameters.h"
#include "tests_contrastgenerator.h"
#include "tests_phantom_input.h"
#include "tests_encoding.h"
#include "tests_mr_acquisition_model.h"

#include "all_simulation_tests.h"



void run_tests_auxiliary_testing_functions( void )
{

	bool tests_successful = true;

	tests_successful *= test_aux_test_funs::test_get_serialized_ismrmrd_header();
	tests_successful *= test_aux_test_funs::test_get_mock_acquisition_vector();

	if ( !tests_successful )
	{
		throw std::runtime_error( "The auxiliary testing functions tests failed.");
	}
	else
	{
		std::cout<< "The auxiliary testing functions tests succeeded." <<std::endl;
	}


}


void run_tests_auxiliary_input_output( void )
{

	bool tests_successful = true;

	test_aux_io::test_write_ndarray_to_raw();

	if ( !tests_successful )
	{
		throw std::runtime_error( "The auxiliary testing functions tests failed.");
	}
	else
	{
		std::cout<< "The auxiliary testing functions tests succeeded." <<std::endl;
	}
}



void run_tests_tissueparameters(void)
{
	bool tests_successful = true;


	// call every test here
	tests_successful *= test_allocate_MRTissueParameter_successful();
	tests_successful *= test_allocate_PETTissueParameter_successful();
	tests_successful *= test_allocate_TissueParameter_successful();
	
	tests_successful *= test_get_MRTissueParameter_from_ptree();
	tests_successful *= test_get_PETTissueParameter_from_ptree();

	//tests_successful *= test_exception_throw_if_node_not_exists();

	tests_successful *= test_read_TissueParameter_label_from_xml(XML_TEST_PATH);
	
	tests_successful *= test_check_label_uniqueness_fails();
	tests_successful *= test_check_label_uniqueness_true();
	

	if ( !tests_successful )
	{
		throw std::runtime_error( "The tissueparameters tests failed.");
	}
	else
	{
		std::cout<< "The tissueparameters tests succeeded." <<std::endl;
	}
}

void run_tests_contrastgenerator(void)
{
	bool tests_successful = true;

	// tlm tests
	tests_successful *= test_tlm::test_get_filepath_tissue_parameter_xml();
	tests_successful *= test_tlm::test_get_labels_array();
	tests_successful *=	test_tlm::test_get_segmentation_dimensions();

	tests_successful *=	test_tlm::test_assign_tissue_parameters_label_found();
	tests_successful *= test_tlm::test_assign_tissue_parameters_label_not_found();

	tests_successful *= test_tlm::test_map_labels_to_tissue_from_xml();


	// contgen tests
	tests_successful *= test_contgen::test_mr_constructor();
	tests_successful *= test_contgen::test_mr_set_rawdata_header();

	tests_successful *=	test_contgen::test_map_flash_contrast();
	tests_successful *=	test_contgen::test_mr_map_contrast_dim_check();

	test_contgen::test_mr_map_contrast_application_to_xcat();

	if ( !tests_successful )
	{
		throw std::runtime_error( "The contrastgenerator tests failed.");
	}
	else
	{
		std::cout<< "The contrastgenerator tests succeeded." <<std::endl;
	}
}


void run_tests_phantom_input( void )
{
	bool tests_successful = true;

	// insert tests
	tests_successful *= test_read_h5_segmentation_correct_dims(H5_PHANTOM_TEST_PATH);
	tests_successful *= test_read_h5_segmentation_correct_content(H5_PHANTOM_TEST_PATH);
	
	test_read_h5_segmentation_for_xcat_input_check(H5_XCAT_PHANTOM_PATH);
		
	if ( !tests_successful )
	{
		throw std::runtime_error( "The h5 file reader tests failed." );
	}
	else
	{
		std::cout<< "The phantom input tests succeeded" << std::endl;
	}


}

void run_tests_encoding( void )
{

	bool tests_successful = true;

	tests_successful *= test_enc::test_cube_input();
	tests_successful *= test_cart_enc::test_sample_fourier_space();
	

	if ( !tests_successful )
	{
		throw std::runtime_error( "The encoding tests failed." );
	}
	else
	{
		std::cout<< "The encoding tests succeeded" << std::endl;
	}


}


void run_tests_mr_acquisition_model( void )
{

	bool tests_successful = true;

	tests_successful *= tests_mracqmod::test_fwd_method();



	if ( !tests_successful )
	{
		throw std::runtime_error( "The acquisition model tests failed." );
	}
	else
	{
		std::cout<< "The acquisition model tests succeeded" << std::endl;
	}


}













