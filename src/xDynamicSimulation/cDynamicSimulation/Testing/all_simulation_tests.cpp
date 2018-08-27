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
#include "tests_dynamics.h"
#include "tests_dynamicsimulation.h"
#include "tests_noisegenerator.h"
#include "tests_dynsim_deformer.h"

#include "all_simulation_tests.h"


void run_tests_auxiliary_testing_functions( void )
{

	bool tests_successful = true;

	tests_successful *= test_aux_test_funs::test_get_serialized_ismrmrd_header();
	tests_successful *= test_aux_test_funs::test_get_mock_acquisition_vector();
	tests_successful *= test_aux_test_funs::test_get_mock_csm();
	tests_successful *= test_aux_test_funs::test_get_mock_coildata_as_cfimage();
	tests_successful *= test_aux_test_funs::test_get_mock_ismrmrd_image_with_cube();
	tests_successful *= test_aux_test_funs::test_get_mock_pet_contrast_generator();



	if ( !tests_successful )
	{
		throw std::runtime_error( "The auxiliary testing functions tests failed.");
	}
	else
	{
		std::cout<< "The auxiliary testing functions tests succeeded." <<std::endl;
	}


}

void run_tests_dynamics( void )
{
	
	bool tests_successful = true;


	tests_successful *= test_dynamic::test_is_in_bin();
	tests_successful *= test_dynamic::test_intersect_mr_acquisition_data();

	tests_successful *= test_dynamic::test_linear_interpolate_signal();
	tests_successful *= test_dynamic::test_get_set_bins();

	// tests_successful *= test_dynamic::test_bin_mr_acquisitions();

	tests_successful *= test_dynamic::test_motion_dynamic_counter();

	// tests_successful *= test_dynamic::test_motion_dynamic_temp_folder_setup();
	tests_successful *= test_dynamic::test_motion_dynamic_set_motion_fields();
	tests_successful *= test_dynamic::test_motion_dynamic_write_motion_fields();


	if ( !tests_successful )
	{
		throw std::runtime_error( "The dynamics tests failed." );
	}
	else
	{
		std::cout<< "The dynamics tests succeeded" << std::endl;
	}	
}


void run_tests_dynamic_simulation( void )
{

	bool tests_successful = true;
	
	// tests_successful *= tests_mr_dynsim::test_acquisitionsvector_memory_management();

	//tests_successful *= test_lin_combi_gen::test_get_all_combinations();

	// tests_successful *= tests_mr_dynsim::test_constructor();
	// tests_mr_dynsim::test_extract_hdr_information();
	tests_successful *= tests_mr_dynsim::test_simulate_contrast_dynamics( );


	// tests_successful *= test_pet_dynsim::test_constructor();
	// tests_successful *= test_pet_dynsim::set_template_acquisition_data();
	// tests_successful *= test_pet_dynsim::test_simulate_dynamics();

	if ( !tests_successful )
	{
		throw std::runtime_error( "The dynamic simulation tests failed." );
	}
	else
	{
		std::cout<< "The dynamic simulation tests succeeded" << std::endl;
	}



}

void run_tests_noise_generator( void )
{

	bool tests_successful = true;

	// tests_successful *= test_noisegen::test_add_poisson_noise();
	tests_successful *= test_noisegen::test_add_gaussian_noise();

	if ( !tests_successful )
	{
		throw std::runtime_error( "The noise generator tests failed." );
	}
	else
	{
		std::cout<< "The noise generator tests succeeded" << std::endl;
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


void run_tests_auxiliary_input_output( void )
{

	bool tests_successful = true;

	// test_aux_io::test_write_ndarray_to_raw();
	test_aux_io::test_write_ismrmrd_image_to_analyze();

	// tests_successful *= test_aux_io::test_read_acquisitions_vector_number_consistency();	


	if ( !tests_successful )
	{
		throw std::runtime_error( "The auxiliary input output functions tests failed.");
	}
	else
	{
		std::cout<< "The auxiliary input output functions tests succeeded." <<std::endl;
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
	
	tests_successful *= test_TissueParameter_algebra();



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

	// // tlm tests
	// tests_successful *= test_tlm::test_get_filepath_tissue_parameter_xml();
	// tests_successful *= test_tlm::test_get_labels_array();
	// tests_successful *=	test_tlm::test_get_segmentation_dimensions();

	// tests_successful *=	test_tlm::test_assign_tissue_parameters_label_found();
	// tests_successful *= test_tlm::test_assign_tissue_parameters_label_not_found();

	// tests_successful *= test_tlm::test_map_labels_to_tissue_from_xml();
	// tests_successful *= test_tlm::test_replace_petmr_tissue_parameters();


	// // mr contgen tests
	// tests_successful *= test_contgen::test_mr_constructor();
	// tests_successful *= test_contgen::test_mr_set_rawdata_header();

	// tests_successful *=	test_contgen::test_map_flash_contrast();
	// tests_successful *=	test_contgen::test_mr_map_contrast_dim_check();


	// test_contgen::test_match_output_dims_to_headerinfo();
	// test_contgen::test_mr_map_contrast_application_to_xcat();
	test_contgen::test_replace_petmr_tissue_parameters_in_xcat();

	// // pet contgen tests
	// tests_successful *=	test_contgen::test_pet_constructor();
	// tests_successful *= test_contgen::test_pet_map_contrast();
	// tests_successful *= test_contgen::test_pet_map_attenuation(); 
	// tests_successful *= test_contgen::test_set_template_image_from_file();

	// test_contgen::test_pet_map_contrast_application_to_xcat();


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

	std::cout<< "outcommented tests failing since file disappeared" << std::endl;
	// tests_successful *= test_read_h5_segmentation_correct_dims(H5_PHANTOM_TEST_PATH);
	// tests_successful *= test_read_h5_segmentation_correct_content(H5_PHANTOM_TEST_PATH);
	
	test_read_h5_segmentation_for_xcat_input_check(H5_XCAT_PHANTOM_PATH);
	tests_successful *= test_read_h5_motionfields();

	
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




void run_tests_dynsim_deformer( void )
{
	// bool tests_successful = true;
	bool tests_successful = true;

	tests_successful *=	DynSimDeformerTester::test_deform_contrast_generator();


	if ( !tests_successful )
	{
		throw std::runtime_error( "The dynsim deformer tests failed." );
	}
	else
	{
		std::cout<< "The dynsim deformer tests succeeded" << std::endl;
	}

}








