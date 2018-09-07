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
		int const N = 1;
		int const M = 2;
		int const L = 3;

		DimensionsType dims;
		dims.push_back(N);
		dims.push_back(M);
		dims.push_back(L);



		LinearCombiGenerator lcg( dims );

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

	MRDynamicSimulation mr_dyn_sim( mr_cont );

	return true;


	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}


}

void tests_mr_dynsim::test_extract_hdr_information( void )
{
	try
	{
	
	MRContrastGenerator mr_cont_gen = aux_test::get_mock_mr_contrast_generator();

	MRDynamicSimulation mr_dyn_sim( mr_cont_gen );

	mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );

	mr_dyn_sim.extract_hdr_information();


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

bool tests_mr_dynsim::test_acquisitionsvector_memory_management( void )
{

	AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( ISMRMRD_H5_TEST_PATH );


	uint const num_reps = 100;
	for(uint i_rep=0; i_rep<num_reps; i_rep++)
	{


		AcquisitionsVector temp_dummy_vector;
		temp_dummy_vector.copy_acquisitions_info(all_acquis);

		for(size_t i_acq=0; i_acq<all_acquis.number(); i_acq++)
		{
			ISMRMRD::Acquisition acq;
			all_acquis.get_acquisition(i_acq, acq);
			temp_dummy_vector.append_acquisition(acq);
		}

		std::cout << "Iteration number: " << i_rep << std::endl;
		std::cout << "Number count in vector: " << temp_dummy_vector.number() << std::endl;
	}

	return true;

}

bool tests_mr_dynsim::test_simulate_contrast_dynamics( void )
{
	try
	{
	
		ISMRMRD::NDArray< unsigned int > segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );
		
		// float const test_noise_width = 0.1;
		// mr_dyn_sim.set_noise_width( test_noise_width );
		float const test_SNR = 15;
		mr_dyn_sim.set_SNR(test_SNR);
		
		int const num_simul_states_first_dyn = 10;
		int const num_simul_states_second_dyn = 10;


		ContrastDynamic first_cont_dyn(num_simul_states_first_dyn), second_cont_dyn(num_simul_states_second_dyn);

		std::vector<LabelType> first_dynamic_labels = {1, 3, 4};	
		for(int i=0; i<first_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << first_dynamic_labels[i] << " to first dynamic." << std::endl;
			first_cont_dyn.add_dynamic_label(first_dynamic_labels[i]);
		}

		std::vector<LabelType> second_dynamic_labels = {5, 6, 7, 8, 36, 37};	
		for(int i=0; i<second_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << second_dynamic_labels[i] << " to second dynamic." << std::endl;
			second_cont_dyn.add_dynamic_label(second_dynamic_labels[i]);
		}


		auto extreme_tissue_params = aux_test::get_mock_contrast_signal_extremes();

		first_cont_dyn.set_parameter_extremes(extreme_tissue_params.first, extreme_tissue_params.second);

		auto second_extremes_0 = extreme_tissue_params.first;
		auto second_extremes_1 = extreme_tissue_params.second;

		second_extremes_0.mr_tissue_.spin_density_percentH2O_ = 95;
		second_extremes_0.mr_tissue_.t1_miliseconds_ = 1000;
		second_extremes_0.mr_tissue_.t2_miliseconds_= 100;
		
		second_extremes_1.mr_tissue_.spin_density_percentH2O_ = 95;
		second_extremes_1.mr_tissue_.t1_miliseconds_ = 500;
		second_extremes_1.mr_tissue_.t2_miliseconds_= 100;

		second_cont_dyn.set_parameter_extremes(second_extremes_0, second_extremes_1);


		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( mr_dyn_sim.get_filename_rawdata() );

		SignalContainer mock_sinus_signal = aux_test::get_mock_sinus_signal(all_acquis);
		SignalContainer mock_ramp_signal = aux_test::get_mock_contrast_signal(all_acquis);


	 	first_cont_dyn.set_dyn_signal( mock_sinus_signal );
	 	second_cont_dyn.set_dyn_signal( mock_ramp_signal );

		first_cont_dyn.bin_mr_acquisitions( all_acquis );
		second_cont_dyn.bin_mr_acquisitions( all_acquis );

		mr_dyn_sim.add_dynamic( first_cont_dyn );
		mr_dyn_sim.add_dynamic( second_cont_dyn );
		
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);
		mr_dyn_sim.simulate_dynamics();

		mr_dyn_sim.write_simulation_results( FILENAME_MR_CONTRAST_DYNSIM );

		return true;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}

}


bool tests_mr_dynsim::test_simulate_motion_dynamics( )
{
	try
	{
		ISMRMRD::NDArray< unsigned int > segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );
		
		// float const test_noise_width = 0.1;
		// mr_dyn_sim.set_noise_width( test_noise_width );
		float const test_SNR = 15;
		mr_dyn_sim.set_SNR(test_SNR);
		
		int const num_simul_cardiac_states = 10;
		int const num_simul_resp_states = 10;
		
		MotionDynamic cardiac_dyn(num_simul_cardiac_states), resp_dyn(num_simul_resp_states);

		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( mr_dyn_sim.get_filename_rawdata() );
		SignalContainer mock_sinus_signal = aux_test::get_mock_sinus_signal(all_acquis);
		SignalContainer mock_ramp_signal = aux_test::get_mock_contrast_signal(all_acquis);

	 	cardiac_dyn.set_dyn_signal( mock_sinus_signal );
	 	cardiac_dyn.bin_mr_acquisitions( all_acquis );

	 	resp_dyn.set_dyn_signal( mock_ramp_signal );
	 	resp_dyn.bin_mr_acquisitions( all_acquis );
		
		// auto motion_fields = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		auto cardiac_motion_fields = read_cardiac_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		auto resp_motion_fields = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		
		cardiac_dyn.set_displacement_fields( cardiac_motion_fields, true );
		resp_dyn.set_displacement_fields( resp_motion_fields, false );


		mr_dyn_sim.add_dynamic( cardiac_dyn );
		mr_dyn_sim.add_dynamic( resp_dyn );
		
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);
		mr_dyn_sim.simulate_dynamics();

		mr_dyn_sim.write_simulation_results( FILENAME_MR_MOTION_DYNSIM );

		return true;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}

}

bool tests_mr_dynsim::test_simulate_simultaneous_motion_contrast_dynamics()
{
	try
	{	
		ISMRMRD::NDArray< unsigned int > segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );
		
		// float const test_noise_width = 0.1;
		// mr_dyn_sim.set_noise_width( test_noise_width );
		float const test_SNR = 15;
		mr_dyn_sim.set_SNR(test_SNR);
		

		int const num_simul_motion_dyn = 5;
		
		MotionDynamic first_motion_dyn(num_simul_motion_dyn), second_motion_dyn( num_simul_motion_dyn );

		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( mr_dyn_sim.get_filename_rawdata() );
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);


		// SignalContainer mock_motion_signal = aux_test::get_mock_sinus_signal(all_acquis);
		SignalContainer mock_motion_signal = aux_test::get_mock_contrast_signal(all_acquis);

		// SETTING UP MOTION DYNAMICS ########################################################################

	 	first_motion_dyn.set_dyn_signal( mock_motion_signal );
	 	first_motion_dyn.bin_mr_acquisitions( all_acquis );
		
		second_motion_dyn.set_dyn_signal( mock_motion_signal );
	 	second_motion_dyn.bin_mr_acquisitions( all_acquis );

		auto motion_fields = read_cardiac_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		first_motion_dyn.set_displacement_fields( motion_fields, true );

		motion_fields = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		second_motion_dyn.set_displacement_fields( motion_fields, false );

		mr_dyn_sim.add_dynamic( first_motion_dyn );
		mr_dyn_sim.add_dynamic( second_motion_dyn );


		// SETTING UP CONRAST DYNAMICS ########################################################################

		int const num_simul_states_first_contrast_dyn = 10;
		int const num_simul_states_second_contrast_dyn = 10;


		ContrastDynamic first_cont_dyn(num_simul_states_first_contrast_dyn), second_cont_dyn(num_simul_states_second_contrast_dyn);

		std::vector<LabelType> first_dynamic_labels = {1, 3, 4};	
		for(int i=0; i<first_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << first_dynamic_labels[i] << " to first dynamic." << std::endl;
			first_cont_dyn.add_dynamic_label(first_dynamic_labels[i]);
		}

		std::vector<LabelType> second_dynamic_labels = {5, 6, 7, 8, 36, 37};	
		for(int i=0; i<second_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << second_dynamic_labels[i] << " to second dynamic." << std::endl;
			second_cont_dyn.add_dynamic_label(second_dynamic_labels[i]);
		}


		auto extreme_tissue_params = aux_test::get_mock_contrast_signal_extremes();

		first_cont_dyn.set_parameter_extremes(extreme_tissue_params.first, extreme_tissue_params.second);

		auto second_extremes_0 = extreme_tissue_params.first;
		auto second_extremes_1 = extreme_tissue_params.second;

		second_extremes_0.mr_tissue_.spin_density_percentH2O_ = 95;
		second_extremes_0.mr_tissue_.t1_miliseconds_ = 1000;
		second_extremes_0.mr_tissue_.t2_miliseconds_= 100;
		
		second_extremes_1.mr_tissue_.spin_density_percentH2O_ = 95;
		second_extremes_1.mr_tissue_.t1_miliseconds_ = 500;
		second_extremes_1.mr_tissue_.t2_miliseconds_= 100;

		second_cont_dyn.set_parameter_extremes(second_extremes_0, second_extremes_1);

		SignalContainer mock_contrast_signal = aux_test::get_mock_contrast_signal(all_acquis);

		first_cont_dyn.set_dyn_signal( mock_contrast_signal );
	 	second_cont_dyn.set_dyn_signal( mock_contrast_signal );

		first_cont_dyn.bin_mr_acquisitions( all_acquis );
		second_cont_dyn.bin_mr_acquisitions( all_acquis );

		mr_dyn_sim.add_dynamic( first_cont_dyn );
		mr_dyn_sim.add_dynamic( second_cont_dyn );
		
		// ####################################################################################################

		clock_t t;
		t = clock();
		mr_dyn_sim.simulate_dynamics();
		t = clock() - t;

		std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;
		mr_dyn_sim.write_simulation_results( FILENAME_MR_MOTION_CONTRAST_DYNSIM );

		return true;
	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}
}



bool test_pet_dynsim::test_constructor()
{

	try
	{
		PETContrastGenerator pet_cont_gen = aux_test::get_mock_pet_contrast_generator();
		PETDynamicSimulation pet_dyn_sim( pet_cont_gen );

		return true;


	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}

}




bool test_pet_dynsim::set_template_acquisition_data()
{

	try
	{
		PETContrastGenerator pet_cont_gen = aux_test::get_mock_pet_contrast_generator();
		PETDynamicSimulation pet_dyn_sim( pet_cont_gen );

		pet_dyn_sim.set_filename_rawdata( PET_TEMPLATE_ACQUISITION_DATA_PATH );

		pet_dyn_sim.set_template_acquisition_data();

		return true;


	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}

}


bool test_pet_dynsim::test_simulate_dynamics()
{

	try
	{
		PETContrastGenerator pet_cont_gen = aux_test::get_mock_pet_contrast_generator();

		PETDynamicSimulation pet_dyn_sim( pet_cont_gen );
		
		pet_dyn_sim.set_filename_rawdata( PET_TEMPLATE_ACQUISITION_DATA_PATH );
		
		pet_dyn_sim.simulate_dynamics();
		
		pet_dyn_sim.write_simulation_results(FILENAME_DYNSIM_PET);

		return true;


	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}

}




