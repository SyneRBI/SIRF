/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


// Files collecting all tests for one specific module.


#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <stdexcept>

#include "test_input_filenames.h"

#include "app_mracquisitiondata.h"
#include "tests_auxiliary_testing_functions.h"
#include "tests_auxiliary_input_output.h"
#include "tests_tissueparameters.h"
#include "tests_contrastgenerator.h"
#include "tests_phantom_input.h"
#include "tests_dynamics.h"
#include "tests_dynamicsimulation.h"
#include "tests_noisegenerator.h"
#include "tests_dynsim_deformer.h"
#include "tests_c_interface.h"


bool run_apps(void)
{
	// std::string const fname_ismrmrd = std::string(SHARED_FOLDER_PATH) + "/PublicationData/FatWaterQuantification/Output/5DMotion/output_grpe_mri_simulation_motion_type_cardiorespiratory__num_motion_states_10_x_10.h5";
	std::string const fname_ismrmrd = std::string(SHARED_FOLDER_PATH) + "/PublicationData/FatWaterQuantification/Output/4DMotion/Cardiac/output_grpe_mri_simulation_motion_type_cardiac_num_motion_states_10.h5";
	apps_johannesmayer::omitt_first_acquisition(fname_ismrmrd);
}


bool run_tests_auxiliary_testing_functions(void )
{

	std::cout<< "Running " << __FUNCTION__ << std::endl;		
	
	try
    {
        bool tests_successful = true;
		int i=0;
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		tests_successful *= test_aux_test_funs::test_get_serialized_ismrmrd_header();
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		tests_successful *= test_aux_test_funs::test_get_mock_acquisition_vector();
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		tests_successful *= test_aux_test_funs::test_get_mock_ismrmrd_image_with_cube();
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		tests_successful *= test_aux_test_funs::test_get_mock_pet_contrast_generator();
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		tests_successful *= test_aux_test_funs::test_get_mock_sawtooth_signal();
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		tests_successful *= test_aux_test_funs::test_get_mock_gaussian_csm();
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;

		return tests_successful;

    }
    catch(std::runtime_error const &e){
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
	catch(const std::exception &error){
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
		throw;
	}
	catch(...){
        std::cerr << "An unknown exception was caught in "<< __FUNCTION__ << std::endl;
		throw;
	}
}

bool run_tests_dynamics( void )
{

	bool tests_successful = true;
	std::vector< bool > dyn_tests;
	std::cout << "start ----------------------------------------------------" <<std::endl;
	dyn_tests.push_back(test_dynamic::test_is_in_bin());

	dyn_tests.push_back(test_surrogateprocessor::test_linear_interpolate_signal());
	dyn_tests.push_back(test_binprocessor::test_get_set_bins());

	dyn_tests.push_back(test_contrastprocessor::test_mr_contrast_motion_dyn_get_num_simul_states());

	dyn_tests.push_back(test_motionprocessor::test_motion_dynamic_counter());
	dyn_tests.push_back(test_motionprocessor::test_motion_dynamic_temp_folder_setup());
	dyn_tests.push_back(test_motionprocessor::test_motion_dynamic_save_gt_deformations());	
	dyn_tests.push_back(test_motionprocessor::test_motion_dynamic_set_motion_fields());	
	dyn_tests.push_back(test_motionprocessor::test_motion_dynamic_prep_motion_fields());	
	dyn_tests.push_back(test_motionprocessor::test_motion_dynamic_temp_interpolate_dvfs());
	dyn_tests.push_back(test_motionprocessor::test_nonisotropic_mvf_resampling () );
	// dyn_tests.push_back(test_dynamic::test_mvf_vs_pet_img_quarternions());
	
	dyn_tests.push_back(test_dynamic::test_intersect_mr_acquisition_data());
	
	dyn_tests.push_back(test_dynamic::test_bin_mr_acquisitions());
	dyn_tests.push_back(test_dynamic::test_bin_pet_time_interval());
	
	std::cout << "dynamics test results = ";
	for( size_t i=0; i<dyn_tests.size(); i++)
	{
		std::cout << dyn_tests[i] << " / ";
		tests_successful *= dyn_tests[i];
	}
	std::cout << std::endl;

	
	if ( !tests_successful )
	{
		std::stringstream ss_msg;
		ss_msg << "Running " << __FUNCTION__ << " failed.";
		throw std::runtime_error( ss_msg.str() );
	}
	else
	{
		std::cout<< "Running " << __FUNCTION__ << " succeeded.";
		return tests_successful;
	}
}


bool run_tests_dynamic_simulation( void )
{
	bool tests_successful = true;
	std::vector< bool > mr_dynsim_tests;

	// mr_dynsim_tests.push_back(test_lin_combi_gen::test_get_all_combinations());
	// mr_dynsim_tests.push_back(tests_datageneration::read_write_h5_filecontent());
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_constructor());
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_simulate_statics());
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_simulate_dynamics());
	mr_dynsim_tests.push_back(tests_mr_dynsim::test_simulate_5d_motion_dynamics());
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_simulate_rpe_acquisition());
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_dce_acquisition());
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_4d_mri_acquisition());
	// mr_dynsim_tests.push_back(tests_mr_dynsim::test_5d_mri_acquisition());

	std::cout << "mr dynamic simulation test results = ";
	for( size_t i=0; i<mr_dynsim_tests.size(); i++)
	{
		std::cout << mr_dynsim_tests[i] << " / ";
		tests_successful *= mr_dynsim_tests[i];
	}
	std::cout << std::endl;

	std::vector< bool > pet_dynsim_tests;

	// pet_dynsim_tests.push_back(test_pet_dynsim::test_constructor());
	// pet_dynsim_tests.push_back(test_pet_dynsim::set_template_acquisition_data());
	// pet_dynsim_tests.push_back(test_pet_dynsim::test_simulate_statics());
	// pet_dynsim_tests.push_back(test_pet_dynsim::test_simulate_motion_dynamics());
	// pet_dynsim_tests.push_back(test_pet_dynsim::test_4d_pet_acquisition());
	// pet_dynsim_tests.push_back(test_pet_dynsim::test_5d_pet_acquisition());
	
	std::cout << "pet dynamic simulation test results = ";
	for( size_t i=0; i<pet_dynsim_tests.size(); i++)
	{
		std::cout << pet_dynsim_tests[i] << " / ";
		tests_successful *= pet_dynsim_tests[i];
	}
	std::cout << std::endl;

	
	if ( !tests_successful )
	{
		std::stringstream ss_msg;
		ss_msg << "Running " << __FUNCTION__ << " failed.";
		throw std::runtime_error( ss_msg.str() );
	}
	else
	{
		std::cout<< "Running " << __FUNCTION__ << " succeeded.";
		return tests_successful;
	}
}


bool run_tests_noise_generator( void )
{
	bool tests_successful = true;
	std::vector<bool> noise_tests;

	noise_tests.push_back(test_noisegen::test_add_poisson_noise());
	noise_tests.push_back(test_noisegen::test_add_gaussian_noise());

	std::cout << "Results " << __FUNCTION__ << " = ";

	for( size_t i=0; i<noise_tests.size(); i++)
	{
		std::cout << noise_tests[i] << " / ";
		tests_successful *= noise_tests[i];
	}

	std::cout << std::endl;

	return tests_successful;	
	
}

bool run_tests_auxiliary_input_output( void )
{
	std::cout<< "Running " << __FUNCTION__ << std::endl;		

	try
    {
        bool tests_successful = true;
		int i=0;
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		test_aux_io::test_write_ndarray_to_raw();
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		test_aux_io::test_write_ismrmrd_image_to_analyze();
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		tests_successful *= test_aux_io::test_read_acquisitions_vector_number_consistency();	
		std::cout << "#:" << ++i << "-------------------------------------" << std::endl;
		tests_successful *= test_aux_io::test_read_single_column_txt();

		return tests_successful;

    }
    catch(std::runtime_error const &e){
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
	catch(const std::exception &error){
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
		throw;
	}
	catch(...){
        std::cerr << "An unknown exception was caught in "<< __FUNCTION__ << std::endl;
		throw;
	}
}

bool run_tests_tissueparameters(void)
{
	std::cout<< "Running " << __FUNCTION__ << std::endl;		

	try
    {
		bool tests_successful = true;
        int i=0;
		// call every test here
		tests_successful *= test_allocate_MRTissueParameter_successful();
		tests_successful *= test_allocate_PETTissueParameter_successful();
		tests_successful *= test_allocate_TissueParameter_successful();
		tests_successful *= test_get_MRTissueParameter_from_ptree();
		tests_successful *= test_get_PETTissueParameter_from_ptree();
		test_exception_throw_if_node_not_exists();
		tests_successful *= test_read_TissueParameter_label_from_xml(XML_TEST_PATH);
		tests_successful *= test_check_label_uniqueness_fails();
		tests_successful *= test_check_label_uniqueness_true();
		tests_successful *= test_TissueParameter_algebra();			

		return tests_successful;

    }
    catch(std::runtime_error const &e){
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
	catch(const std::exception &error){
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
		throw;
	}
	catch(...){
        std::cerr << "An unknown exception was caught in "<< __FUNCTION__ << std::endl;
		throw;
	}
}

bool run_tests_contrastgenerator(void)
{
try{
		
	int i=0;
	bool tests_successful = true;
	std::vector< bool > tlm_tests, abstract_contgen_tests, mr_contgen_tests, pet_contgen_tests;
	
	// tlm tests

	tlm_tests.push_back( test_tlm::test_get_filepath_tissue_parameter_xml() );
	tlm_tests.push_back( test_tlm::test_get_labels_array() );
	tlm_tests.push_back( test_tlm::test_get_segmentation_dimensions() );
	tlm_tests.push_back( test_tlm::test_assign_tissue_parameters_label_found() );
	tlm_tests.push_back( test_tlm::test_assign_tissue_parameters_label_not_found() );
	tlm_tests.push_back( test_tlm::test_map_labels_to_tissue_from_xml() );
	tlm_tests.push_back( test_tlm::test_replace_petmr_tissue_parameters() );

	std::cout << "#### #### #### Tissue-Label-Mapper test results = ";
	for( size_t i=0; i<tlm_tests.size(); i++)
	{
		std::cout << tlm_tests[i] << " / ";
		tests_successful *= tlm_tests[i];
	}
	std::cout << std::endl;
	// abstract contgent tests

	abstract_contgen_tests.push_back( test_contgen::test_get_tissue_parameter() );

	std::cout << "#### #### #### Abstract contrast generator test results = ";
	for( size_t i=0; i<abstract_contgen_tests.size(); i++)
	{
		std::cout << abstract_contgen_tests[i] << " / ";
		tests_successful *= abstract_contgen_tests[i];
	}
	std::cout << std::endl;

	// mr contgen tests
	mr_contgen_tests.push_back( test_contgen::test_mr_map_contrast_dim_check() );
	mr_contgen_tests.push_back( test_contgen::test_mr_constructor() );
	mr_contgen_tests.push_back( test_contgen::test_mr_set_rawdata_header() );
	mr_contgen_tests.push_back( test_contgen::test_map_flash_contrast() );
	
	test_contgen::test_mr_map_external_contrast_to_xcat();
	test_contgen::test_mr_map_contrast_application_to_xcat();
	test_contgen::test_replace_petmr_tissue_parameters_in_xcat();
	test_contgen::test_get_signal_for_tissuelabel_in_xcat();

	std::cout << "#### #### #### MR constrast generator test results = ";
	
	for( size_t i=0; i<mr_contgen_tests.size(); i++)
	{
		std::cout << mr_contgen_tests[i] << " / ";
		tests_successful *= mr_contgen_tests[i];
	}
	std::cout << std::endl;

	// pet contgen tests

	pet_contgen_tests.push_back( test_contgen::test_pet_constructor() );
	pet_contgen_tests.push_back( test_contgen::test_pet_map_contrast() );
	pet_contgen_tests.push_back( test_contgen::test_pet_map_attenuation() ); 
	pet_contgen_tests.push_back( test_contgen::test_set_template_image_from_file() );
	pet_contgen_tests.push_back( test_contgen::test_resample_to_template_image() );
	test_contgen::test_pet_map_contrast_application_to_xcat();

	std::cout << "#### #### #### PET Contrast generator test results = ";
	for( size_t i=0; i<pet_contgen_tests.size(); i++)
	{
		std::cout << pet_contgen_tests[i] << " / ";
		tests_successful *= pet_contgen_tests[i];
	}
	std::cout << std::endl;

	return tests_successful;

	}
    catch(std::runtime_error const &e){
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
	catch(const std::exception &error){
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
		throw;
	}
	catch(...){
        std::cerr << "An unknown exception was caught in "<< __FUNCTION__ << std::endl;
		throw;
	}
}


bool run_tests_phantom_input( void )
{
	std::cout<< "Running " << __FUNCTION__ << std::endl;		

	try
    {

		bool tests_successful = true;
		std::vector< bool > test_results{};

		test_results.push_back(test_read_1D_dataset_from_h5(H5_XCAT_PHANTOM_PATH));
		test_results.push_back(test_read_geometrical_info_from_h5( H5_XCAT_PHANTOM_PATH ));
		test_results.push_back(test_read_segmentation_to_nifti( H5_XCAT_PHANTOM_PATH ));
		test_results.push_back(test_read_motionfield_to_nifti(  H5_XCAT_PHANTOM_PATH ));

		std::cout << "#### #### #### " << __FUNCTION__ << " test results = ";
		for( size_t i=0; i<test_results.size(); i++)
		{
			std::cout << test_results[i] << " / ";
			tests_successful *= test_results[i];
		}
		std::cout << std::endl;

		return tests_successful;

    }
    catch(std::runtime_error const &e){
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
	catch(const std::exception &error){
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
		throw;
	}
	catch(...){
        std::cerr << "An unknown exception was caught in "<< __FUNCTION__ << std::endl;
		throw;
	}
}


bool run_tests_dynsim_deformer( void )
{
	std::cout<< "Running " << __FUNCTION__ << std::endl;		
	try
    {
		bool tests_successful = true;
		std::vector< bool > test_results{};

		test_results.push_back(DynSimDeformerTester::test_nifti_data_deformation());
		test_results.push_back(DynSimDeformerTester::test_deform_mr_contrast_generator());
		test_results.push_back(DynSimDeformerTester::test_deform_pet_contrast_generator());
		test_results.push_back(DynSimDeformerTester::test_motion_of_MotionDynamics());


		std::cout << "#### #### #### " << __FUNCTION__ << " test results = ";
		for( size_t i=0; i<test_results.size(); i++)
		{
			std::cout << test_results[i] << " / ";
			tests_successful *= test_results[i];
		}
		
		std::cout << std::endl;

		return tests_successful;

    }
    catch(std::runtime_error const &e){
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
	catch(const std::exception &error){
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
		throw;
	}
	catch(...){
        std::cerr << "An unknown exception was caught in "<< __FUNCTION__ << std::endl;
		throw;
	}
}


bool run_tests_c_interface( void )
{
	std::cout<< "Running " << __FUNCTION__ << std::endl;		
	try
    {
		bool tests_successful = true;
		std::vector< bool > test_results{};

		test_results.push_back(test_simulation_interface::test_bin_data_from_handle());
	
		std::cout << "#### #### #### " << __FUNCTION__ << " test results = ";
		for( size_t i=0; i<test_results.size(); i++)
		{
			std::cout << test_results[i] << " / ";
			tests_successful *= test_results[i];
		}
		
		std::cout << std::endl;

		return tests_successful;

    }
    catch(std::runtime_error const &e){
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
	catch(const std::exception &error){
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
		throw;
	}
	catch(...){
        std::cerr << "An unknown exception was caught in "<< __FUNCTION__ << std::endl;
		throw;
	}
}


int main ( int argc, char* argv[])
{

	try{
		std::cout << "Starting Simulation C++ tests... " <<std::endl;
		bool ok = true;

		if(argc > 1)
			fprintf(stdout, "Please do not pass any arguments. This just runs test code.");

		// ok *= run_tests_auxiliary_testing_functions();
		// ok *= run_tests_auxiliary_input_output();
		// ok *= run_tests_tissueparameters();
		// ok *= run_tests_contrastgenerator();
		// ok *= run_tests_phantom_input();
		// ok *= run_tests_noise_generator();
		// ok *= run_tests_dynamics();
		// ok *= run_tests_c_interface();
		// ok *= run_tests_dynsim_deformer();
		ok *= run_tests_dynamic_simulation();
				
		if(ok)
			return EXIT_SUCCESS;	
		else
			return EXIT_FAILURE;
	}
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
	
}


