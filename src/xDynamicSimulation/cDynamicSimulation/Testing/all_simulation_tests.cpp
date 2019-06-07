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
#include "tests_volume_orientator.h"
#include "tests_memory_usage.h" 

#include "all_simulation_tests.h"


void run_tests_auxiliary_testing_functions( void )
{

	bool tests_successful = true;

	// tests_successful *= test_aux_test_funs::test_get_serialized_ismrmrd_header();
	// tests_successful *= test_aux_test_funs::test_get_mock_acquisition_vector();
	// tests_successful *= test_aux_test_funs::test_get_mock_csm();
	// tests_successful *= test_aux_test_funs::test_get_mock_coildata_as_cfimage();
	// tests_successful *= test_aux_test_funs::test_get_mock_ismrmrd_image_with_cube();
	// tests_successful *= test_aux_test_funs::test_get_mock_pet_contrast_generator();
	// tests_successful *= test_aux_test_funs::test_get_mock_sawtooth_signal();
	tests_successful *= test_aux_test_funs::test_get_mock_gaussian_csm();


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
	std::vector< bool > dyn_tests;
	std::cout << "start ----------------------------------------------------" <<std::endl;
	// std::cout << "1 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_is_in_bin());

	// std::cout << "2 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_intersect_mr_acquisition_data());

	// std::cout << "3 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_linear_interpolate_signal());

	// std::cout << "4 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_get_set_bins());

	// std::cout << "5 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_bin_mr_acquisitions());

	// std::cout << "6 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_motion_dynamic_counter());

	// std::cout << "7 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_motion_dynamic_temp_folder_setup());

	// std::cout << "7.1 --------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_motion_dynamic_save_gt_deformations());	

	// std::cout << "8 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_motion_dynamic_set_motion_fields());	

	// std::cout << "10 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_motion_dynamic_prep_motion_fields());	

	// std::cout << "11 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_motion_dynamic_temp_interpolate_dvfs());

	// std::cout << "12 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_mvf_vs_pet_img_quarternions());

	// std::cout << "13 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_mr_contrast_motion_dyn_get_num_simul_states());
	
	// std::cout << "14 ----------------------------------------------------" <<std::endl;
	dyn_tests.push_back(test_dynamic::test_bin_pet_time_interval());

	std::cout << "15 ----------------------------------------------------" <<std::endl;
	// dyn_tests.push_back(test_dynamic::test_nonisotropic_mvf_resampling () );
	
	std::cout << "end ----------------------------------------------------" <<std::endl;


	std::cout << "dynamics test results = ";
	for( size_t i=0; i<dyn_tests.size(); i++)
	{
		std::cout << dyn_tests[i] << " / ";
		tests_successful *= dyn_tests[i];
	}
	std::cout << std::endl;

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
	std::vector< bool > mr_dynsim_tests;

	std::cout << "start ----------------------------------------------------" <<std::endl;
	std::cout << "MR 1 ----------------------------------------------------" <<std::endl;
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_acquisitionsvector_memory_management());

	std::cout << "MR 2 ----------------------------------------------------" <<std::endl;
	// mr_dynsim_tests.push_back(test_lin_combi_gen::test_get_all_combinations());

	std::cout << "MR 3 ----------------------------------------------------" <<std::endl;
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_constructor());

	std::cout << "MR void test 1 ----------------------------------------------------" <<std::endl;
	// tests_mr_dynsim::test_extract_hdr_information();

	std::cout << "MR 4 ----------------------------------------------------" <<std::endl;
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_simulate_dynamics());

	std::cout << "MR 5 ----------------------------------------------------" <<std::endl;
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_simulate_rpe_acquisition());

	std::cout << "MR 6 ----------------------------------------------------" <<std::endl;
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_dce_acquisition());

	std::cout << "MR 7 ----------------------------------------------------" <<std::endl;
	mr_dynsim_tests.push_back(tests_mr_dynsim::test_4d_mri_acquisition());

	
	std::cout << "mr dynamic simulation test results = ";
	for( size_t i=0; i<mr_dynsim_tests.size(); i++)
	{
		std::cout << mr_dynsim_tests[i] << " / ";
		tests_successful *= mr_dynsim_tests[i];
	}
	std::cout << std::endl;


	std::vector< bool > pet_dynsim_tests;

	std::cout << "PET 1 ----------------------------------------------------" <<std::endl;
	// pet_dynsim_tests.push_back(test_pet_dynsim::test_constructor());

	std::cout << "PET 2 ----------------------------------------------------" <<std::endl;
	// pet_dynsim_tests.push_back(test_pet_dynsim::set_template_acquisition_data());

	std::cout << "PET 3 ----------------------------------------------------" <<std::endl;
	// pet_dynsim_tests.push_back(test_pet_dynsim::test_simulate_statics());

	std::cout << "PET 4 ----------------------------------------------------" <<std::endl;
	// pet_dynsim_tests.push_back(test_pet_dynsim::test_simulate_motion_dynamics());

	std::cout << "PET 5 ----------------------------------------------------" <<std::endl;
	// pet_dynsim_tests.push_back(test_pet_dynsim::test_4d_pet_acquisition());


	std::cout << "pet dynamic simulation test results = ";
	for( size_t i=0; i<pet_dynsim_tests.size(); i++)
	{
		std::cout << pet_dynsim_tests[i] << " / ";
		tests_successful *= pet_dynsim_tests[i];
	}
	std::cout << std::endl;


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
	std::cout << "Starting aux io tests" << std::endl;
	bool tests_successful = true;

	// test_aux_io::test_write_ndarray_to_raw();
	// test_aux_io::test_write_ismrmrd_image_to_analyze();

	// tests_successful *= test_aux_io::test_read_acquisitions_vector_number_consistency();	
	tests_successful *= test_aux_io::test_read_single_column_txt();

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
	std::vector< bool > tlm_tests, abstract_contgen_tests, mr_contgen_tests, pet_contgen_tests;

	
	// // tlm tests
	// std::cout << "----------------------------------------------------" <<std::endl;
	// tlm_tests.push_back( test_tlm::test_get_filepath_tissue_parameter_xml() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// tlm_tests.push_back( test_tlm::test_get_labels_array() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// tlm_tests.push_back( test_tlm::test_get_segmentation_dimensions() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// tlm_tests.push_back( test_tlm::test_assign_tissue_parameters_label_found() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// tlm_tests.push_back( test_tlm::test_assign_tissue_parameters_label_not_found() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// tlm_tests.push_back( test_tlm::test_map_labels_to_tissue_from_xml() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// tlm_tests.push_back( test_tlm::test_replace_petmr_tissue_parameters() );
	// std::cout << "----------------------------------------------------" <<std::endl;

	// std::cout << "tlm test results = ";
	// for( size_t i=0; i<tlm_tests.size(); i++)
	// {
	// 	std::cout << tlm_tests[i] << " / ";
	// 	tests_successful *= tlm_tests[i];
	// }
	// std::cout << std::endl;

	// abstract contgent tests
	// std::cout << "----------------------------------------------------" <<std::endl;
	// abstract_contgen_tests.push_back( test_contgen::test_get_tissue_parameter() );
		for( size_t i=0; i<abstract_contgen_tests.size(); i++)
	{
		std::cout << abstract_contgen_tests[i] << " / ";
		tests_successful *= abstract_contgen_tests[i];
	}
	std::cout << std::endl;

	// // mr contgen tests
	// std::cout << "----------------------------------------------------" <<std::endl;
	// mr_contgen_tests.push_back( test_contgen::test_mr_constructor() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// mr_contgen_tests.push_back( test_contgen::test_mr_set_rawdata_header() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// mr_contgen_tests.push_back( test_contgen::test_map_flash_contrast() );
	// std::cout << "----------------------------------------------------" <<std::endl;
	// mr_contgen_tests.push_back( test_contgen::test_mr_map_contrast_dim_check() );
	// std::cout << "----------------------------------------------------" <<std::endl;

	// test_contgen::test_match_output_dims_to_headerinfo();
	// std::cout << "----------------------------------------------------" <<std::endl;
	// test_contgen::test_mr_map_contrast_application_to_xcat();
	// std::cout << "----------------------------------------------------" <<std::endl;
	// test_contgen::test_replace_petmr_tissue_parameters_in_xcat();
	// std::cout << "----------------------------------------------------" <<std::endl;
	// test_contgen::test_get_signal_for_tissuelabel_in_xcat();
	// std::cout << "----------------------------------------------------" <<std::endl;
	// std::cout << "mr contgen test results = ";
	// for( size_t i=0; i<mr_contgen_tests.size(); i++)
	// {
	// 	std::cout << mr_contgen_tests[i] << " / ";
	// 	tests_successful *= mr_contgen_tests[i];
	// }
	// std::cout << std::endl;

	// pet contgen tests
	// std::cout << "----------------------------------------------------" <<std::endl;
	// pet_contgen_tests.push_back( test_contgen::test_pet_constructor() );
	std::cout << "----------------------------------------------------" <<std::endl;
	// pet_contgen_tests.push_back( test_contgen::test_pet_map_contrast() );
	std::cout << "----------------------------------------------------" <<std::endl;
	// pet_contgen_tests.push_back( test_contgen::test_pet_map_attenuation() ); 
	std::cout << "----------------------------------------------------" <<std::endl;
	// pet_contgen_tests.push_back( test_contgen::test_set_template_image_from_file() );
	std::cout << "----------------------------------------------------" <<std::endl;
	pet_contgen_tests.push_back( test_contgen::test_resample_to_template_image() );
	
	std::cout << "----------------------------------------------------" <<std::endl;
	// test_contgen::test_pet_map_contrast_application_to_xcat();
	std::cout << "----------------------------------------------------" <<std::endl;

	// std::cout << "pet contgen test results = ";
	// for( size_t i=0; i<pet_contgen_tests.size(); i++)
	// {
	// 	std::cout << pet_contgen_tests[i] << " / ";
	// 	tests_successful *= pet_contgen_tests[i];
	// }
	std::cout << std::endl;

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


	// test_read_1D_dataset_from_h5(H5_PHANTOM_TEST_PATH);
	// test_read_geometrical_info_from_h5( H5_PHANTOM_TEST_PATH );
	// test_read_segmentation_to_nifti( H5_PHANTOM_TEST_PATH );
	test_read_motionfield_to_nifti(  H5_XCAT_PHANTOM_PATH );

	// tests_successful *= test_read_h5_segmentation_correct_dims(H5_XCAT_PHANTOM_PATH);
	// tests_successful *= test_read_h5_segmentation_correct_content(H5_XCAT_PHANTOM_PATH);
	
	// test_read_h5_segmentation_for_xcat_input_check(H5_XCAT_PHANTOM_PATH);
	// tests_successful *= test_read_h5_motionfields();

	
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
	// tests_successful *= CartesianEncodingTester::test_sample_fourier_space();
	

	// tests_successful *= RPETrajectoryPreparationTester::test_get_set_trajectory();
	// tests_successful *= RPETrajectoryPreparationTester::test_get_result_container();

	tests_successful *= RPETester::test_sample_fourier_space();

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
	bool tests_successful = true;

	std::cout << " Start -------------------------- " <<std::endl;

	std::cout << " 0-------------------------- " <<std::endl;
	// tests_successful *= DynSimDeformerTester::test_nifti_data_deformation();

	// std::cout << " 1-------------------------- " <<std::endl;
	tests_successful *=	DynSimDeformerTester::test_deform_contrast_generator();

	std::cout << " 2 -------------------------- " <<std::endl;
	// tests_successful *= DynSimDeformerTester::test_SIRFImageDataDeformation_memory_behavior();

	std::cout << " 3 -------------------------- " <<std::endl;
	// tests_successful *= DynSimDeformerTester::test_deform_pet_contrast_generator();

	std::cout << " 4-------------------------- " <<std::endl;
	// tests_successful *= DynSimDeformerTester::test_motion_of_MotionDynamics();

	std::cout << " End -------------------------- " <<std::endl;

	if ( !tests_successful )
	{
		throw std::runtime_error( "The dynsim deformer tests failed." );
	}
	else
	{
		std::cout<< "The dynsim deformer tests succeeded" << std::endl;
	}

}




void run_tests_volume_orientator( void )
{
	bool tests_successful = true;

	tests_successful *= sirf::aVolumeOrientatorTester::test_reorient_image();

	if ( !tests_successful )
	{
		throw std::runtime_error( "The volume orientator tests failed." );
	}
	else
	{
		std::cout<< "The volume orientator tests succeeded" << std::endl;
	}

}


void run_tests_memory_usage( void )
{
	bool tests_successful = true;

	tests_memory::test_acquisition_memory();
	// tests_memory::test_downsizing_acquisition_memory();
	// tests_successful *= tests_memory::test_acquisition_vector_memory();
	// tests_successful *=tests_memory::test_acquisition_vector_ordering_memory();
	// tests_successful *= tests_memory::test_ndarray_memory_managment();
	// tests_successful *= tests_memory::tests_resizing_acquisition_memory();
	// tests_successful *= tests_memory::tests_VD_h5_file_content();
	if ( !tests_successful )
	{
		throw std::runtime_error( "The memory usage tests failed." );
	}
	else
	{
		std::cout<< "The memory usage tests succeeded" << std::endl;
	}



}