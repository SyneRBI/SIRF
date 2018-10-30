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

#include "gadgetron_data_containers.h"

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
		
		// size_t const NRad = 64;
		// size_t const NAng = 64;
		RPETrajectoryContainer rpe_traj;
		auto sptr_traj = std::make_shared< RPETrajectoryContainer >( rpe_traj );
		mr_dyn_sim.set_trajectory( sptr_traj );
		
		float const test_SNR = 15;
		mr_dyn_sim.set_SNR(test_SNR);
		
		int const num_simul_states_myocardium_dyn = 10;
		int const num_simul_states_blood_dyn = 1;



		MRContrastDynamic myocardium_cont_dyn(num_simul_states_myocardium_dyn), blood_cont_dyn(num_simul_states_blood_dyn);

		std::vector<LabelType> myocardium_dynamic_labels = {1, 3, 4};	
		for(int i=0; i<myocardium_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << myocardium_dynamic_labels[i] << " to myocardium dynamic." << std::endl;
			myocardium_cont_dyn.add_dynamic_label(myocardium_dynamic_labels[i]);
		}

		std::vector<LabelType> blood_dynamic_labels = {5, 6, 7, 8, 36, 37};	
		for(int i=0; i<blood_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << blood_dynamic_labels[i] << " to vascular dynamic." << std::endl;
			blood_cont_dyn.add_dynamic_label(blood_dynamic_labels[i]);
		}


		auto extreme_tissue_params = aux_test::get_mock_contrast_signal_extremes();

		myocardium_cont_dyn.set_parameter_extremes(extreme_tissue_params.first, extreme_tissue_params.second);

		auto blood_extremes_0 = extreme_tissue_params.first;
		auto blood_extremes_1 = extreme_tissue_params.second;

		blood_extremes_0.mr_tissue_.spin_density_percentH2O_ = 95;
		blood_extremes_0.mr_tissue_.t1_miliseconds_ = 1000;
		blood_extremes_0.mr_tissue_.t2_miliseconds_= 100;
		
		blood_extremes_1.mr_tissue_.spin_density_percentH2O_ = 95;
		blood_extremes_1.mr_tissue_.t1_miliseconds_ = 500;
		blood_extremes_1.mr_tissue_.t2_miliseconds_= 100;

		blood_cont_dyn.set_parameter_extremes(blood_extremes_0, blood_extremes_1);


		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( mr_dyn_sim.get_filename_rawdata() );

		SignalContainer mock_ramp_signal = aux_test::get_generic_contrast_inflow_signal(all_acquis);

	 	myocardium_cont_dyn.set_dyn_signal( mock_ramp_signal );
	 	blood_cont_dyn.set_dyn_signal( mock_ramp_signal );

		myocardium_cont_dyn.bin_mr_acquisitions( all_acquis );
		blood_cont_dyn.bin_mr_acquisitions( all_acquis );


		mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (myocardium_cont_dyn) );
		// mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (blood_cont_dyn) );
		
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);

		clock_t t;
		t = clock();
		mr_dyn_sim.simulate_dynamics();
		t = clock() - t;

		std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;
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
	

		auto data_dims = segmentation_labels.getDims();
		
		std::vector< size_t > vol_dims{data_dims[0], data_dims[1], data_dims[2]}; 
		
		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( csm );

		std::string const traj_name = "Cartesian";

		if( traj_name == "ITLGCRPE") 
		{
			RPEInterleavedGoldenCutTrajectoryContainer rpe_traj;
			auto sptr_traj = std::make_shared< RPEInterleavedGoldenCutTrajectoryContainer >( rpe_traj );
			mr_dyn_sim.set_trajectory( sptr_traj );
		}
		else if (traj_name == "RPE")
		{
			RPETrajectoryContainer sf_traj;
			auto sptr_traj = std::make_shared< RPETrajectoryContainer >( sf_traj );
			mr_dyn_sim.set_trajectory( sptr_traj );
		}
		else if (traj_name == "SF")
		{
			SFTrajectoryContainer sf_traj;
			auto sptr_traj = std::make_shared< SFTrajectoryContainer >( sf_traj );
			mr_dyn_sim.set_trajectory( sptr_traj );
		}
		else if (traj_name == "GCSF")
		{
			GoldenAngleSFTrajectoryContainer gcsf_traj;
			auto sptr_traj = std::make_shared< GoldenAngleSFTrajectoryContainer >( gcsf_traj );
			mr_dyn_sim.set_trajectory( sptr_traj );
		}


		float const test_SNR = 150;
		mr_dyn_sim.set_SNR(test_SNR);

		int const num_simul_cardiac_states = 10;
		int const num_simul_resp_states = 10;
		
		MRMotionDynamic cardiac_dyn(num_simul_cardiac_states), resp_dyn(num_simul_resp_states);

		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( mr_dyn_sim.get_filename_rawdata() );

		SignalContainer mock_cardiac_signal = aux_test::get_generic_contrast_inflow_signal(all_acquis);
		SignalContainer mock_respiratory_signal = aux_test::get_generic_contrast_inflow_signal(all_acquis);

	 	cardiac_dyn.set_dyn_signal( mock_cardiac_signal );
	 	cardiac_dyn.bin_mr_acquisitions( all_acquis );

	 	resp_dyn.set_dyn_signal( mock_respiratory_signal );
	 	resp_dyn.bin_mr_acquisitions( all_acquis );
		
		auto cardiac_motion_fields = read_cardiac_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		auto resp_motion_fields = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		
		cardiac_dyn.set_displacement_fields( cardiac_motion_fields, true );
		resp_dyn.set_displacement_fields( resp_motion_fields, false );

		// mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> (cardiac_dyn) );
		mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> (resp_dyn) );
		
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);
		
		clock_t t;
		t = clock();
		mr_dyn_sim.simulate_dynamics();
		t = clock() - t;

		std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;

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
		
		
		float const test_SNR = 150;
		mr_dyn_sim.set_SNR(test_SNR);
		
		int const num_simul_motion_dyn = 10;

		// RPETrajectoryContainer rpe_traj;
		// auto sptr_traj = std::make_shared< RPETrajectoryContainer >( rpe_traj );
		// mr_dyn_sim.set_trajectory( sptr_traj );
		

		MRMotionDynamic cardiac_motion_dyn(num_simul_motion_dyn), respiratory_motion_dyn( num_simul_motion_dyn );

		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( mr_dyn_sim.get_filename_rawdata() );
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);


		SignalContainer mock_cardiac_signal = aux_test::get_generic_cardiac_signal(all_acquis);
		SignalContainer mock_respiratory_signal = aux_test::get_generic_respiratory_signal(all_acquis);
		

		// SETTING UP MOTION DYNAMICS ########################################################################

	 	cardiac_motion_dyn.set_dyn_signal( mock_cardiac_signal );
	 	cardiac_motion_dyn.bin_mr_acquisitions( all_acquis );
		
		respiratory_motion_dyn.set_dyn_signal( mock_respiratory_signal );
	 	respiratory_motion_dyn.bin_mr_acquisitions( all_acquis );

		auto motion_fields = read_cardiac_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		cardiac_motion_dyn.set_displacement_fields( motion_fields, true );

		motion_fields = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		respiratory_motion_dyn.set_displacement_fields( motion_fields, false );


		mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> (cardiac_motion_dyn ));
		mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> (respiratory_motion_dyn ));


		// SETTING UP CONRAST DYNAMICS ########################################################################

		int const num_simul_states_myocardium_contrast_dyn = 10;
		int const num_simul_states_blood_contrast_dyn = 10;



		MRContrastDynamic myocardium_cont_dyn(num_simul_states_myocardium_contrast_dyn), blood_cont_dyn(num_simul_states_blood_contrast_dyn);

		std::vector<LabelType> myocardium_dynamic_labels = {1, 3, 4};	
		for(int i=0; i<myocardium_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << myocardium_dynamic_labels[i] << " to myocardium dynamic." << std::endl;
			myocardium_cont_dyn.add_dynamic_label(myocardium_dynamic_labels[i]);
		}

		std::vector<LabelType> blood_dynamic_labels = {5, 6, 7, 8, 36, 37};	
		for(int i=0; i<blood_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << blood_dynamic_labels[i] << " to vascular dynamic." << std::endl;
			blood_cont_dyn.add_dynamic_label(blood_dynamic_labels[i]);
		}


		auto extreme_tissue_params = aux_test::get_mock_contrast_signal_extremes();

		myocardium_cont_dyn.set_parameter_extremes(extreme_tissue_params.first, extreme_tissue_params.second);

		auto blood_extremes_0 = extreme_tissue_params.first;
		auto blood_extremes_1 = extreme_tissue_params.second;

		blood_extremes_0.mr_tissue_.spin_density_percentH2O_ = 95;
		blood_extremes_0.mr_tissue_.t1_miliseconds_ = 1000;
		blood_extremes_0.mr_tissue_.t2_miliseconds_= 100;
		
		blood_extremes_1.mr_tissue_.spin_density_percentH2O_ = 95;
		blood_extremes_1.mr_tissue_.t1_miliseconds_ = 500;
		blood_extremes_1.mr_tissue_.t2_miliseconds_= 100;

		blood_cont_dyn.set_parameter_extremes(blood_extremes_0, blood_extremes_1);

		SignalContainer mock_contrast_signal = aux_test::get_generic_contrast_inflow_signal(all_acquis);

		myocardium_cont_dyn.set_dyn_signal( mock_contrast_signal );
	 	blood_cont_dyn.set_dyn_signal( mock_contrast_signal );

		myocardium_cont_dyn.bin_mr_acquisitions( all_acquis );
		blood_cont_dyn.bin_mr_acquisitions( all_acquis );

		mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (myocardium_cont_dyn) );
		mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (blood_cont_dyn) );
		
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

bool tests_mr_dynsim::test_simulate_rpe_acquisition()
{
	
	using sirf::RPETrajectoryContainer;


	try
	{	
		ISMRMRD::NDArray< unsigned int > segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );
	
		data_io::write_raw( std::string(SHARED_FOLDER_PATH) + "seg_in_sim", segmentation_labels.begin(), segmentation_labels.getNumberOfElements() );

		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );

		 

		RPEInterleavedGoldenCutTrajectoryContainer rpe_traj;
		auto sptr_traj = std::make_shared< RPEInterleavedGoldenCutTrajectoryContainer >( rpe_traj );
		mr_dyn_sim.set_trajectory( sptr_traj );

		auto data_dims = segmentation_labels.getDims();
		std::vector< size_t > vol_dims{data_dims[0], data_dims[1], data_dims[2]}; 
		
		std::cout << epiph( data_dims[0] ) <<std::endl;
		std::cout << epiph( data_dims[1] ) <<std::endl;
		std::cout << epiph( data_dims[2] ) <<std::endl;

		size_t num_coils = 20;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( csm );

		float const test_SNR = 15;
		mr_dyn_sim.set_SNR(test_SNR);

		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( mr_dyn_sim.get_filename_rawdata() );
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);

		clock_t t;
		t = clock();
		mr_dyn_sim.simulate_statics();
		t = clock() - t;

		std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;
		mr_dyn_sim.write_simulation_results( FILENAME_MR_RPE_SIM );


		return true;
	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}

}


// PET #################################################################################



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


bool test_pet_dynsim::test_simulate_statics()
{

	try
	{
		PETContrastGenerator pet_cont_gen = aux_test::get_mock_pet_contrast_generator();

		PETDynamicSimulation pet_dyn_sim( pet_cont_gen );
		
		pet_dyn_sim.set_filename_rawdata( PET_TEMPLATE_ACQUISITION_DATA_PATH );
		pet_dyn_sim.set_template_image_data( PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH );

		clock_t t;
		t = clock();
		pet_dyn_sim.simulate_statics();
		t = clock() - t;

		std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;

		
		pet_dyn_sim.write_simulation_results(FILENAME_STATICSIM_PET);

		return true;


	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}

}



bool test_pet_dynsim::test_simulate_motion_dynamics()
{

	try
	{
		PETContrastGenerator pet_cont_gen = aux_test::get_mock_pet_contrast_generator();

		PETDynamicSimulation pet_dyn_sim( pet_cont_gen );

		pet_dyn_sim.set_output_filename_prefix("/media/sf_SharedFolder/CCPPETMR/ISMRMSim/Output/DynamicPET/4DResp/pet_dyn_simul");
		
		pet_dyn_sim.set_filename_rawdata( PET_TEMPLATE_ACQUISITION_DATA_PATH );
		pet_dyn_sim.set_template_image_data( PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH );
		
		int const num_simul_cardiac_states = 10;
		PETMotionDynamic  resp_dyn(num_simul_cardiac_states);


		SignalContainer resp_sig = data_io::read_surrogate_signal( std::string(TIME_POINTS_PATH), std::string(RESP_SIGNAL_PATH));
		
		auto first_pt = resp_sig[0];
		auto last_pt = resp_sig[ resp_sig.size()-1 ];

		float min_time_ms = first_pt.first;
		float tot_time_ms = last_pt.first - first_pt.first;

		std::cout << "total time ms: " << tot_time_ms <<std::endl;

		for( size_t i=0; i<resp_sig.size(); i++)
		{
			auto curr_sig_pt = resp_sig[i];	
			curr_sig_pt.first = 25 * (curr_sig_pt.first - min_time_ms)/tot_time_ms;
			resp_sig[i] = curr_sig_pt;
		}

	 	first_pt = resp_sig[0];
		last_pt = resp_sig[ resp_sig.size()-1 ];

		min_time_ms = first_pt.first;
		tot_time_ms = last_pt.first - first_pt.first;

		std::cout << "total time ms: " << tot_time_ms <<std::endl;
	 	resp_dyn.set_dyn_signal( resp_sig );

	 	resp_dyn.bin_total_time_interval( TimeBin(min_time_ms,tot_time_ms) );
		TimeBin total_time(0, 25);

		resp_dyn.bin_total_time_interval( total_time );
		

		auto resp_motion_fields = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		resp_dyn.set_displacement_fields( resp_motion_fields, false );
		
		pet_dyn_sim.add_dynamic( std::make_shared<PETMotionDynamic> (resp_dyn) );
		pet_dyn_sim.simulate_dynamics( tot_time_ms );

		return true;


	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}

}




