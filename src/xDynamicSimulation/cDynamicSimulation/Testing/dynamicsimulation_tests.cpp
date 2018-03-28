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

	tests_successful *= test_tlm::test_set_get_filepath_tissue_parameter_xml();
	tests_successful *= test_tlm::test_set_get_labels_array();

	tests_successful *= test_tlm::test_tlm_constructor();

	tests_successful *=	test_tlm::test_assign_tissue_parameters_label_found();
	tests_successful *= test_tlm::test_assign_tissue_parameters_label_not_found();

	


	if ( !tests_successful )
	{
		throw std::runtime_error( "The tissueparameters tests failed.");
	}
	else
	{
		std::cout<< "The tissueparameters tests succeeded." <<std::endl;
	}
}


void run_tests_phantom_input( void )
{
	bool tests_successful = true;

	// insert tests
	tests_successful *= test_read_h5_segmentation_correct_dims(H5_TEST_PATH);
	tests_successful *= test_read_h5_segmentation_correct_content(H5_TEST_PATH);
	
	
	if ( !tests_successful )
	{
		throw std::runtime_error( "The h5 file reader tests failed." );
	}
	else
	{
		std::cout<< "The h5 file reader tests succeeded" << std::endl;
	}


}