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

#include "sirf/cGadgetron/gadgetron_data_containers.h"

#include "sirf/cDynamicSimulation/dynamicsimulation_x.h"
#include "auxiliary_testing_functions.h"
#include "sirf/cDynamicSimulation/phantom_input.h"

#include "tests_dynamicsimulation.h"


using namespace std;
using namespace sirf;
using namespace ISMRMRD;


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


bool tests_mr_dynsim::test_simulate_dynamics()
{
	std::cout << "Running function " <<__FUNCTION__ <<" .!" <<std::endl;

	try
	{	
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );
		
		auto data_dims = segmentation_labels.get_dimensions();
		
		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
		
		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( csm );


		// std::string const traj_name = "ITLGCRPE";
		std::string const traj_name = "Cartesian";

		if( traj_name == "ITLGCRPE") 
		{
			RPEInterleavedGoldenCutTrajectoryContainer rpe_traj;
			auto sptr_traj = std::make_shared< RPEInterleavedGoldenCutTrajectoryContainer >( rpe_traj );
			mr_dyn_sim.set_trajectory( sptr_traj );
		}


		AcquisitionsVector all_acquis;
		all_acquis.read( mr_dyn_sim.get_filename_rawdata() );
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);

		float const test_SNR = 15;
		size_t const noise_label = 13;
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );
		
		int const num_simul_motion_dyn = 5;
		
		MRMotionDynamic cardiac_motion_dyn(num_simul_motion_dyn), respiratory_motion_dyn( num_simul_motion_dyn );

		SignalContainer mock_cardiac_signal = aux_test::get_mock_sawtooth_signal(all_acquis, 1000);
		SignalContainer mock_respiratory_signal = aux_test::get_mock_sinus_signal(all_acquis, 3000);
		
		// SETTING UP MOTION DYNAMICS ########################################################################

	 // 	cardiac_motion_dyn.set_dyn_signal( mock_cardiac_signal );
	 // 	cardiac_motion_dyn.bin_mr_acquisitions( all_acquis );
		
		// auto cardiac_motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		// cardiac_motion_dyn.set_displacement_fields( cardiac_motion_fields, true );

		// respiratory_motion_dyn.set_dyn_signal( mock_respiratory_signal );
	 // 	respiratory_motion_dyn.bin_mr_acquisitions( all_acquis );

		// auto resp_motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		// respiratory_motion_dyn.set_displacement_fields( resp_motion_fields, false );


		// mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> (cardiac_motion_dyn ));
		// mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> (respiratory_motion_dyn ));


		// SETTING UP CONRAST DYNAMICS ########################################################################

		// int const num_simul_states_myocardium_contrast_dyn = 2;
		// int const num_simul_states_blood_contrast_dyn = 2;

		// MRContrastDynamic myocardium_cont_dyn(num_simul_states_myocardium_contrast_dyn), blood_cont_dyn(num_simul_states_blood_contrast_dyn);

		// std::vector<LabelType> myocardium_dynamic_labels = {1, 2, 3, 4};	
		// for(int i=0; i<myocardium_dynamic_labels.size(); i++)
		// {
		// 	std::cout << "Adding label " << myocardium_dynamic_labels[i] << " to myocardium dynamic." << std::endl;
		// 	myocardium_cont_dyn.add_dynamic_label(myocardium_dynamic_labels[i]);
		// }

		// std::vector<LabelType> blood_dynamic_labels = {5, 6, 7, 8, 36, 37};	
		// for(int i=0; i<blood_dynamic_labels.size(); i++)
		// {
		// 	std::cout << "Adding label " << blood_dynamic_labels[i] << " to vascular dynamic." << std::endl;
		// 	blood_cont_dyn.add_dynamic_label(blood_dynamic_labels[i]);
		// }


		// auto extreme_tissue_params = aux_test::get_mock_contrast_signal_extremes();

		// myocardium_cont_dyn.set_parameter_extremes(extreme_tissue_params.first, extreme_tissue_params.second);

		// auto blood_extremes_0 = extreme_tissue_params.first;
		// auto blood_extremes_1 = extreme_tissue_params.second;

		// blood_extremes_0.mr_tissue_.spin_density_percentH2O_ = 95;
		// blood_extremes_0.mr_tissue_.t1_miliseconds_ = 1000;
		// blood_extremes_0.mr_tissue_.t2_miliseconds_= 100;
		
		// blood_extremes_1.mr_tissue_.spin_density_percentH2O_ = 95;
		// blood_extremes_1.mr_tissue_.t1_miliseconds_ = 500;
		// blood_extremes_1.mr_tissue_.t2_miliseconds_= 100;

		// blood_cont_dyn.set_parameter_extremes(blood_extremes_0, blood_extremes_1);

		// SignalContainer myocard_contrast_signal = aux_test::get_generic_contrast_inflow_signal(all_acquis);
		// SignalContainer blood_contrast_signal = aux_test::get_generic_contrast_in_and_outflow_signal( all_acquis );

		// myocardium_cont_dyn.set_dyn_signal( myocard_contrast_signal );
	 // 	myocardium_cont_dyn.bin_mr_acquisitions( all_acquis );

	 // 	blood_cont_dyn.set_dyn_signal( blood_contrast_signal );
		// blood_cont_dyn.bin_mr_acquisitions( all_acquis );

		// mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (myocardium_cont_dyn) );
		// mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (blood_cont_dyn) );
		
		// ####################################################################################################

		clock_t t;
		t = clock();
		mr_dyn_sim.simulate_dynamics();
		t = clock() - t;

		std::cout << "Storing ground truth motion information" << std::endl;
		mr_dyn_sim.save_ground_truth_displacements();



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

		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );
		
		auto data_dims = segmentation_labels.get_dimensions();


		RPEInterleavedGoldenCutTrajectoryContainer rpe_traj;
		auto sptr_traj = std::make_shared< RPEInterleavedGoldenCutTrajectoryContainer >( rpe_traj );
		mr_dyn_sim.set_trajectory( sptr_traj );

		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
		
		std::cout << epiph( data_dims[0] ) <<std::endl;
		std::cout << epiph( data_dims[1] ) <<std::endl;
		std::cout << epiph( data_dims[2] ) <<std::endl;
		std::cout << epiph( data_dims[3] ) <<std::endl;


		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( csm );

		float const test_SNR = 15;
		size_t const noise_label = 3;
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );

		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( mr_dyn_sim.get_filename_rawdata() );
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);

		clock_t t;
		t = clock();
		mr_dyn_sim.simulate_dynamics();
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



bool tests_mr_dynsim::test_4d_mri_acquisition( void )
{
try
	{	
		bool const do_cardiac_sim = true;
		bool const simulate_data = true;
		bool const store_gt_mvfs = false;

		std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Input/";
		std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Output/MRI/";


		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( std::string(SHARED_FOLDER_PATH) + "/XCATSegmentations/xcat208Cube/xcat_phantom_incl_geomertry_208.h5" );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		mr_dyn_sim.set_filename_rawdata( input_path + "/MRI/meas_MID00241_FID69145_Tho_T1_fast_ismrmrd.h5");
		
		auto data_dims = segmentation_labels.get_dimensions();
		
		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
		
		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( csm );


		RPEInterleavedGoldenCutTrajectoryContainer rpe_traj;
		auto sptr_traj = std::make_shared< RPEInterleavedGoldenCutTrajectoryContainer >( rpe_traj );
		mr_dyn_sim.set_trajectory( sptr_traj );

		AcquisitionsVector all_acquis;
		all_acquis.read( mr_dyn_sim.get_filename_rawdata(), false );
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);

		float const test_SNR = 250;
		size_t const noise_label = 13;
		
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );
		
		int const num_simul_motion_dyn = 24;
		

		// SETTING UP MOTION DYNAMICS ########################################################################

		if( num_simul_motion_dyn > 1)
		{
			MRMotionDynamic motion_dyn( num_simul_motion_dyn );

			std::string const signal_path = input_path + "/SurrogateSignals/";

			std::string fname_timepts, fname_signalpts;

			if( do_cardiac_sim )
			{
				fname_timepts = signal_path + "card_time";
				fname_signalpts = signal_path + "card_signal";
			}
			else
			{
				fname_timepts  = signal_path + "resp_time";
				fname_signalpts = signal_path + "resp_signal";
			}

			SignalContainer motion_signal = data_io::read_surrogate_signal(fname_timepts, fname_signalpts);

		 	motion_dyn.set_dyn_signal( motion_signal );
		 	motion_dyn.bin_mr_acquisitions( all_acquis );

		 	auto binned_acq = motion_dyn.get_binned_mr_acquisitions();
		 	for( int i=0; i<binned_acq.size(); ++i)
		 		std::cout<< "In this bin " << i << " we have " << binned_acq[i].number() << " acquisitions." << std::endl;;

		 	if( do_cardiac_sim )
			{
				auto cardiac_motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
				motion_dyn.set_displacement_fields( cardiac_motion_fields, true );
			}
	 		else
	 		{
				auto resp_motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
				motion_dyn.set_displacement_fields( resp_motion_fields, false );
	 		}

			mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> ( motion_dyn ));
		}
		// ####################################################################################################

		if( simulate_data )
		{
			clock_t t;
			t = clock();
			mr_dyn_sim.simulate_dynamics();
			t = clock() - t;

			std::cout << " TIME FOR 4D MRI SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;

			std::string const motion_type = do_cardiac_sim ? "cardiac" : "respiratory";			

			std::stringstream outname_stream;
			outname_stream << "output_grpe_mri_simulation_" << "motion_type_" << motion_type << "_num_motion_states_" << num_simul_motion_dyn;
			
			std::string const filename_mri_output = output_path + outname_stream.str() + ".h5";

			mr_dyn_sim.write_simulation_results( filename_mri_output );

		}
		if( store_gt_mvfs )
		{
			std::cout << "Storing ground truth motion information" << std::endl;
			mr_dyn_sim.save_ground_truth_displacements();
		}

     	return true;

	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}
}


bool tests_mr_dynsim::test_dce_acquisition( void )
{
	try
	{	
		std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Input/DCE/";
		std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Output/DCE/";


		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( input_path + "Phantoms/xcat_phantom_incl_geomertry_192_dce.h5" );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		mr_dyn_sim.set_filename_rawdata( ISMRMRD_H5_TEST_PATH );
		
		auto data_dims = segmentation_labels.get_dimensions();
		
		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
		
		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( csm );


		RPEInterleavedGoldenCutTrajectoryContainer rpe_traj;
		auto sptr_traj = std::make_shared< RPEInterleavedGoldenCutTrajectoryContainer >( rpe_traj );
		mr_dyn_sim.set_trajectory( sptr_traj );

		AcquisitionsVector all_acquis;
		all_acquis.read( mr_dyn_sim.get_filename_rawdata(), false );
		mr_dyn_sim.set_all_source_acquisitions(all_acquis);

		float const test_SNR = 19;
		size_t const noise_label = 13;
		
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );
		
		int const num_simul_motion_dyn = 1;
		

		// SETTING UP MOTION DYNAMICS ########################################################################
		if( num_simul_motion_dyn > 0)
		{
			MRMotionDynamic respiratory_motion_dyn( num_simul_motion_dyn );

			std::string const filename_resp_timepoints = input_path + "/timepoints_dce_resp_signal";
			std::string const filename_resp_signal = input_path + "/dce_resp_signal";

			SignalContainer respiratory_signal = data_io::read_surrogate_signal(filename_resp_timepoints, filename_resp_signal);

		 	respiratory_motion_dyn.set_dyn_signal( respiratory_signal );
		 	respiratory_motion_dyn.bin_mr_acquisitions( all_acquis );

			auto resp_motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			respiratory_motion_dyn.set_displacement_fields( resp_motion_fields, false );

			mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> ( respiratory_motion_dyn ));
		}

		// SETTING UP CONRAST DYNAMICS ########################################################################

		int const num_contrast_states = 2;
		
		MRContrastDynamic aif_contrast(num_contrast_states), healthy_tissue_contrast(num_contrast_states), lesion_contrast(num_contrast_states);

		std::vector<LabelType> aif_dynamic_labels = {5, 6, 7, 8, 36, 37};	
		for(int i=0; i<aif_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << aif_dynamic_labels[i] << " to AIF dynamic." << std::endl;
			aif_contrast.add_dynamic_label(aif_dynamic_labels[i]);
		}

		std::vector<LabelType> healthy_tissue_dynamic_labels = {13};	
		for(int i=0; i<healthy_tissue_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << healthy_tissue_dynamic_labels[i] << " to healthy tissue dynamic." << std::endl;
			healthy_tissue_contrast.add_dynamic_label(healthy_tissue_dynamic_labels[i]);

		}

		std::vector<LabelType> lesion_dynamic_labels = {74};
		for(int i=0; i<lesion_dynamic_labels.size(); i++)
		{
			std::cout << "Adding label " << lesion_dynamic_labels[i] << " to lesion dynamic." << std::endl;
			lesion_contrast.add_dynamic_label(lesion_dynamic_labels[i]);
		}

		// set up contrast dynamics for vessels
		TissueParameter aif_0 = mr_cont_gen.get_petmr_tissue_parameter( aif_dynamic_labels[0] );
		TissueParameter aif_1 = aif_0;

		aif_0.mr_tissue_.t1_miliseconds_ = 234;

		aif_contrast.set_parameter_extremes( aif_0, aif_1 );

		// set up contrast dynamics for healthy liver
		TissueParameter healthy_tissue_0 = mr_cont_gen.get_petmr_tissue_parameter( healthy_tissue_dynamic_labels[0] );
		TissueParameter healthy_tissue_1 = healthy_tissue_0;

		healthy_tissue_0.mr_tissue_.t1_miliseconds_ = 457;

		healthy_tissue_contrast.set_parameter_extremes( healthy_tissue_0, healthy_tissue_1 );

		// set up contrast dynamics for lesion
		TissueParameter lesion_tissue_0 = mr_cont_gen.get_petmr_tissue_parameter( lesion_dynamic_labels[0] );
		TissueParameter lesion_tissue_1 = lesion_tissue_0;

		lesion_tissue_0.mr_tissue_.t1_miliseconds_ = 432;

		lesion_contrast.set_parameter_extremes( lesion_tissue_0, lesion_tissue_1 );


		std::string const filename_contrast_timepoints = input_path + "/timepoints_dce_contrast_signal";
		
		std::string const filename_aif_t1_ =input_path + "/aif_signal_T1_0_0.23482_T1_1_2.0853";
		std::string const filename_healthy_tissue_t1 = input_path + "/liver_signal_T1_0_0.45683_T1_1_0.80503";
 		std::string const filename_lesion_t1_ = input_path + "/lesion_signal_T1_0_0.43241_T1_1_0.82076";

		SignalContainer aif_dyn_signal = data_io::read_surrogate_signal(filename_contrast_timepoints, filename_aif_t1_);
		SignalContainer healthy_dyn_tissue_signal = data_io::read_surrogate_signal(filename_contrast_timepoints, filename_healthy_tissue_t1);
		SignalContainer lesion_dyn_signal = data_io::read_surrogate_signal(filename_contrast_timepoints, filename_lesion_t1_);


		aif_contrast.set_dyn_signal( aif_dyn_signal );
		healthy_tissue_contrast.set_dyn_signal( healthy_dyn_tissue_signal );
		lesion_contrast.set_dyn_signal( lesion_dyn_signal );


	 	aif_contrast.bin_mr_acquisitions( all_acquis );
		healthy_tissue_contrast.bin_mr_acquisitions( all_acquis );
		lesion_contrast.bin_mr_acquisitions( all_acquis );

		mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (aif_contrast) );
		mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (healthy_tissue_contrast) );
		mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (lesion_contrast) );

		
		// ####################################################################################################

		clock_t t;
		t = clock();
		mr_dyn_sim.simulate_dynamics();
		t = clock() - t;

		std::cout << "Storing ground truth motion information" << std::endl;
		mr_dyn_sim.save_ground_truth_displacements();



		std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;

		std::stringstream outname_stream;
		outname_stream << "output_grpe_dce_simulation_num_motion_states_" << num_simul_motion_dyn <<"_num_contrast_states_" << num_contrast_states;

		std::string const filename_dce_output = output_path + outname_stream.str() + ".h5";

		mr_dyn_sim.write_simulation_results( filename_dce_output );

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
		pet_dyn_sim.set_output_filename_prefix("/media/sf_SharedFolder/CCPPETMR/");
		
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

		// pet_dyn_sim.set_output_filename_prefix("/media/sf_SharedFolder/CCPPETMR/PublicationData/Output/PET/Cardiac4D/pet_dyn_4D_cardiac_simul");
		pet_dyn_sim.set_output_filename_prefix("pet_dyn_4D_cardiac_simul");
		
		pet_dyn_sim.set_filename_rawdata( PET_TEMPLATE_ACQUISITION_DATA_PATH );
		pet_dyn_sim.set_template_image_data( PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH );
		
		// int const num_simul_resp_states = 10;
		int const num_simul_card_states = 24;

		// PETMotionDynamic  resp_dyn(num_simul_resp_states);
		PETMotionDynamic  card_dyn(num_simul_card_states);


		// SignalContainer resp_sig = data_io::read_surrogate_signal( std::string(TIME_POINTS_RESP_PATH), std::string(RESP_SIGNAL_PATH));
			
		// auto first_resp_pt = resp_sig[0];
		// auto last_resp_pt = resp_sig[ resp_sig.size()-1 ];

		// float min_time_resp_ms = first_resp_pt.first;
		// float tot_time_resp_ms = last_resp_pt.first - first_resp_pt.first;

		// for( size_t i=0; i<resp_sig.size(); i++)
		// {
		// 	auto curr_sig_pt = resp_sig[i];	
		// 	curr_sig_pt.first = 25000 * (curr_sig_pt.first - min_time_resp_ms)/tot_time_resp_ms;
		// 	resp_sig[i] = curr_sig_pt;
		// }

		SignalContainer card_sig = data_io::read_surrogate_signal( std::string(TIME_POINTS_CARDIAC_PATH), std::string(CARDIAC_SIGNAL_PATH));

		auto first_card_pt = card_sig[0];
		auto last_card_pt = card_sig[ card_sig.size()-1 ];

		float min_time_card_ms = first_card_pt.first;
		float tot_time_card_ms = last_card_pt.first - first_card_pt.first;


		for( size_t i=0; i<card_sig.size(); i++)
		{
			auto curr_sig_pt = card_sig[i];	
			curr_sig_pt.first = curr_sig_pt.first - min_time_card_ms;
			card_sig[i] = curr_sig_pt;
		}

	 	// resp_dyn.set_dyn_signal( resp_sig );
	 	card_dyn.set_dyn_signal( card_sig );

	 	TimeBin total_time(0, tot_time_card_ms);

	 	// resp_dyn.bin_total_time_interval( total_time );
	 	card_dyn.bin_total_time_interval( total_time );
		

		// auto resp_motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		auto card_motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		
		// resp_dyn.set_displacement_fields( resp_motion_fields, false );
		card_dyn.set_displacement_fields( card_motion_fields, true );
		
		// pet_dyn_sim.add_dynamic( std::make_shared<PETMotionDynamic> (resp_dyn) );
		pet_dyn_sim.add_dynamic( std::make_shared<PETMotionDynamic> (card_dyn) );
		

		pet_dyn_sim.simulate_dynamics( tot_time_card_ms );

		return true;


	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}

}


bool test_pet_dynsim::test_4d_pet_acquisition()
{

	try
	{

		bool const do_cardiac_sim = true;
		bool const simulate_data = false;
		bool const store_gt_mvfs = true;

		std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Input/";
		// std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Output/PET/";

		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		PETContrastGenerator pet_cont_gen( segmentation_labels, XML_XCAT_PATH);

		pet_cont_gen.set_template_image_from_file( PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH );

		PETDynamicSimulation pet_dyn_sim( pet_cont_gen );

		// fix name of output
		stringstream prefix_stream;
		prefix_stream << "pet_dyn_4D_";
		if(do_cardiac_sim){prefix_stream << "cardiac";}else{ prefix_stream << "resp"; } prefix_stream << "_simul";

		pet_dyn_sim.set_output_filename_prefix(prefix_stream.str());
		
		pet_dyn_sim.set_filename_rawdata( PET_TEMPLATE_ACQUISITION_DATA_PATH );
		pet_dyn_sim.set_template_image_data( PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH );
		
		int const num_sim_motion_states = 24;

		float tot_time_ms = 30 * 60 * 1000; // 20 Minute Exam

		if( num_sim_motion_states > 1 )

		{
			PETMotionDynamic motion_dyn( num_sim_motion_states );

			SignalContainer motion_signal;

			stringstream path_time_pts, path_sig_pts;

			path_time_pts  << input_path << "/SurrogateSignals/";
			path_sig_pts << input_path << "/SurrogateSignals/";

			if( do_cardiac_sim )
			{
				path_time_pts << "card_time";
				path_sig_pts  << "card_signal";
			}
			else
			{
				path_time_pts << "resp_time";
				path_sig_pts  << "resp_signal";
			}

			motion_signal = data_io::read_surrogate_signal( path_time_pts.str(), path_sig_pts.str());

			auto first_card_pt = motion_signal[0];
			auto last_card_pt = motion_signal[ motion_signal.size()-1 ];

			float min_time_card_ms = first_card_pt.first;
			tot_time_ms = last_card_pt.first - first_card_pt.first;


			for( size_t i=0; i<motion_signal.size(); i++)
			{
				auto curr_sig_pt = motion_signal[i];	
				curr_sig_pt.first = curr_sig_pt.first - min_time_card_ms;
				motion_signal[i] = curr_sig_pt;
			}

		 	motion_dyn.set_dyn_signal( motion_signal );

		 	TimeBin total_time(0, tot_time_ms);

		 	motion_dyn.bin_total_time_interval( total_time );
			
		 	if( do_cardiac_sim )
			{
				auto motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
				motion_dyn.set_displacement_fields( motion_fields, true );
			}
			else
			{ 		
				auto motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
				motion_dyn.set_displacement_fields( motion_fields, false );
			}
			
			pet_dyn_sim.add_dynamic( std::make_shared<PETMotionDynamic> (motion_dyn) );

		}
		if( simulate_data )
		{
			std::cout << "Simulating Data" << std::endl;
			pet_dyn_sim.simulate_dynamics( tot_time_ms );
			std::cout << "Finished Simulating Data" << std::endl;
		}

		if( store_gt_mvfs )
		{
			std::cout << "Storing GT MVFs" << std::endl;
			pet_dyn_sim.save_ground_truth_displacements();
			std::cout << "Finished Storing GT MVFs" << std::endl;
		}

		return true;

	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}

}