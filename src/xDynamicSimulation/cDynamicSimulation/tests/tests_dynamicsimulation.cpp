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

#include "sirf/Gadgetron/gadgetron_data_containers.h"

#include "sirf/cDynamicSimulation/dynamicsimulation_x.h"
#include "auxiliary_testing_functions.h"
#include "sirf/cDynamicSimulation/phantom_input.h"

#include "tests_dynamicsimulation.h"


using namespace std;
using namespace sirf;
using namespace ISMRMRD;


bool test_lin_combi_gen::test_get_all_combinations( void )
{
	std::cout << " --- Running " << __FUNCTION__ << std::endl;

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
	std::cout << " --- Running " << __FUNCTION__ << std::endl;
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

bool tests_mr_dynsim::test_simulate_dynamics()
{

	std::cout << " --- Running function " <<__FUNCTION__ <<" .!" <<std::endl;

	try
	{	
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);
		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );

		AcquisitionsVector all_acquis;
		all_acquis.read( ISMRMRD_H5_TEST_PATH );
		mr_dyn_sim.set_template_acquisition_data(all_acquis);
		
		auto data_dims = segmentation_labels.get_dimensions();
		
		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
		
		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( std::make_shared<CoilSensitivitiesVector>(csm));

		float const test_SNR = 15;
		size_t const noise_label = 13;
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );

		clock_t t;
		t = clock();
		mr_dyn_sim.simulate_dynamics();
		t = clock() - t;

		std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;

		std::stringstream ss_output_name;
		ss_output_name << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_test_" << __FUNCTION__ << ".h5";
		mr_dyn_sim.write_simulation_results( ss_output_name.str() );

		return true;
	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " << __FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}
}

bool tests_mr_dynsim::test_simulate_rpe_acquisition()
{
	std::cout << " --- Running function " << __FUNCTION__ <<" .!" <<std::endl;
	try
	{	
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		
		auto data_dims = segmentation_labels.get_dimensions();


		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
		
		std::cout << epiph( data_dims[0] ) <<std::endl;
		std::cout << epiph( data_dims[1] ) <<std::endl;
		std::cout << epiph( data_dims[2] ) <<std::endl;
		std::cout << epiph( data_dims[3] ) <<std::endl;


		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( std::make_shared<CoilSensitivitiesVector>(csm));

		float const test_SNR = 15;
		size_t const noise_label = 13;
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );

		AcquisitionsVector all_acquis(ISMRMRD_H5_TEST_PATH);
		mr_dyn_sim.set_template_acquisition_data(all_acquis);

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

bool tests_mr_dynsim::test_5d_mri_acquisition( void )
{
	std::cout << " --- Running function " <<__FUNCTION__ <<" .!" <<std::endl;
	try
		{	
		bool const simulate_data = false;
		bool const store_gt_mvfs = true;

		int const num_simul_motion_states = 10;

		float const test_SNR = 10;
		size_t const noise_label = 13;

		// std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Input/";
		// std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Output/MRI/5DMotion/";

		std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/FatWaterQuantification/Input/";
		std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/FatWaterQuantification/Output/5DMotion/";

		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		std::string fname_rawdata = input_path + "/MR/meas_MID00443_FID81493_3DFatWater_Rpe_Sfl_bSSFP_5min_ismrmrd.h5";
		AcquisitionsVector all_acquis;
		all_acquis.read(fname_rawdata);

		mr_dyn_sim.set_template_acquisition_data(all_acquis);
		auto data_dims = segmentation_labels.get_dimensions();
			
		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
			
		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( std::make_shared<CoilSensitivitiesVector>(csm));

		

					
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );
			
		// SETTING UP MOTION DYNAMICS ########################################################################

		if( num_simul_motion_states > 1)
		{
			MRMotionDynamic card_dyn(num_simul_motion_states), resp_dyn(num_simul_motion_states);

			std::string const signal_path = input_path + "/SurrogateSignals/";

			std::string fname_timepts, fname_signalpts;

			// add card motion
			fname_timepts = signal_path + "card_time";
			fname_signalpts = signal_path + "card_signal";
			SignalContainer card_signal = data_io::read_surrogate_signal(fname_timepts, fname_signalpts);
			card_dyn.set_dyn_signal( card_signal );

			// add resp motion
			fname_timepts  = signal_path + "resp_time";
			fname_signalpts = signal_path + "resp_signal";
			SignalContainer resp_signal = data_io::read_surrogate_signal(fname_timepts, fname_signalpts);
			resp_dyn.set_dyn_signal( resp_signal );

			card_dyn.set_ground_truth_folder_name( output_path + "ground_truth_motionfields_card");
			resp_dyn.set_ground_truth_folder_name( output_path + "ground_truth_motionfields_resp");

			card_dyn.bin_mr_acquisitions( all_acquis );
			resp_dyn.bin_mr_acquisitions( all_acquis );

			auto binned_resp_acq = resp_dyn.get_binned_mr_acquisitions();
			auto binned_card_acq = card_dyn.get_binned_mr_acquisitions();
			for( int i=0; i<binned_resp_acq.size(); ++i)
				std::cout<< "In this resp bin " << i << " we have " << binned_resp_acq[i].number() << " acquisitions." << std::endl;
			for( int i=0; i<binned_resp_acq.size(); ++i)
				std::cout<< "In this card bin " << i << " we have " << binned_card_acq[i].number() << " acquisitions." << std::endl;;
			
			auto cardiac_motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			card_dyn.set_displacement_fields( cardiac_motion_fields, true );
			
			auto resp_motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			resp_dyn.set_displacement_fields( resp_motion_fields, false );

			mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> ( card_dyn ));
			mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> ( resp_dyn ));
		}
		// ####################################################################################################



		if( simulate_data )
		{
			clock_t t;
			t = clock();
			mr_dyn_sim.simulate_dynamics();
			t = clock() - t;

			std::cout << " TIME FOR 5D MRI SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;
			
			std::stringstream outname_stream;
			outname_stream << "output_grpe_mri_simulation_" << "motion_type_cardiorespiratory_"<< "_num_motion_states_" << num_simul_motion_states << "_x_"<< num_simul_motion_states;
			
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


bool tests_mr_dynsim::test_4d_mri_acquisition( void )
{
	try
	{	
		bool const do_cardiac_sim = true;
		bool const simulate_data = true;
		bool const store_gt_mvfs = false;

		int const num_simul_motion_dyn = 4;

		float const test_SNR = 18;
		size_t const noise_label = 13;

		std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Input/";
        std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Output/ComputTimeMeasurement/";

		// std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Output/MRI/5DMotion/";

		// std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/FatWaterQuantification/Input/";
		// std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/FatWaterQuantification/Output/4DMotion/Cardiac/";

		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		std::string fname_rawdata = input_path + "/MRI/meas_MID00241_FID69145_Tho_T1_fast_ismrmrd.h5"; // PETMR

		AcquisitionsVector all_acquis;
		all_acquis.read(fname_rawdata);
		mr_dyn_sim.set_template_acquisition_data(all_acquis);			


		std::vector<float> roi_labels{1,2,3,4,50,72,73};
		std::string const output_prefix_roi = output_path;
		
		auto data_dims = segmentation_labels.get_dimensions();
		
		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
		
		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( std::make_shared<CoilSensitivitiesVector>(csm));



				
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );
		
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

			std::string motion_type_suffix = do_cardiac_sim? "card" : "resp";
			motion_dyn.set_ground_truth_folder_name( output_path + "ground_truth_motionfields_" + motion_type_suffix);


			SignalContainer motion_signal = data_io::read_surrogate_signal(fname_timepts, fname_signalpts);

			// std::cout << "WARNING: CONSTANT SIGNAL: " << std::endl;
			// for( size_t j=0; j< motion_signal.size(); ++j)
			// 	motion_signal[j].second = 0.99;

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


		aux_test::store_roi( segmentation_labels, roi_labels, output_prefix_roi);


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
		bool const simulate_data = false;
		bool const store_gt_mvfs = true;

		std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Input/DCE/";
		std::string const output_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Output/DCE/";

		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( input_path + "Phantoms/xcat_phantom_incl_geomertry_192_dce.h5" );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		MRDynamicSimulation mr_dyn_sim( mr_cont_gen );
		
		AcquisitionsVector all_acquis;
		all_acquis.read(ISMRMRD_H5_TEST_PATH);
		mr_dyn_sim.set_template_acquisition_data(all_acquis);		
		
		auto data_dims = segmentation_labels.get_dimensions();
		std::vector< size_t > vol_dims{(size_t)data_dims[1], (size_t)data_dims[2], (size_t)data_dims[3]}; 
		
		size_t num_coils = 4;
		auto csm = aux_test::get_mock_gaussian_csm(vol_dims, num_coils);
		mr_dyn_sim.set_coilmaps( std::make_shared<CoilSensitivitiesVector>(csm));

		float const test_SNR = 19;
		size_t const noise_label = 13;
		
		mr_dyn_sim.set_SNR(test_SNR);
		mr_dyn_sim.set_noise_label( noise_label );
		
		int const num_simul_motion_dyn = 16;
		
		// SETTING UP MOTION DYNAMICS ########################################################################
		if( num_simul_motion_dyn > 0)
		{
			MRMotionDynamic respiratory_motion_dyn( num_simul_motion_dyn );

			std::string const filename_resp_timepoints = input_path + "/timepoints_dce_resp_signal";
			std::string const filename_resp_signal = input_path + "/dce_resp_signal";

			SignalContainer respiratory_signal = data_io::read_surrogate_signal(filename_resp_timepoints, filename_resp_signal);

		 	respiratory_motion_dyn.set_dyn_signal( respiratory_signal );
		 	respiratory_motion_dyn.bin_mr_acquisitions( all_acquis );

			respiratory_motion_dyn.set_ground_truth_folder_name( output_path + "ground_truth_motionfields");


			auto resp_motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			respiratory_motion_dyn.set_displacement_fields( resp_motion_fields, false );

			mr_dyn_sim.add_dynamic( std::make_shared<MRMotionDynamic> ( respiratory_motion_dyn ));
		}

		// SETTING UP CONRAST DYNAMICS ########################################################################

		int const num_contrast_states = 48;
		
		int const t1_aif_0	= 235;
		int const t1_lesion_0	= 527;
		int const t1_healthy_tissue_0	= 488;


		std::string const filename_contrast_timepoints = input_path + "/timepoints_dce_contrast_signal";
		
		std::string const filename_aif_t1_ =input_path + "/aif_signal_T1_0_0.23482_T1_1_2.0853";
		std::string const filename_healthy_tissue_t1 = input_path + "/liver_signal_T1_0_0.48778_T1_1_0.90637";
 		std::string const filename_lesion_t1_ = input_path + "/lesion_signal_T1_0_0.52744_T1_1_1.2473";

		// ####################################################################################################

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

		aif_0.mr_tissue_.t1_miliseconds_ = t1_aif_0;

		aif_contrast.set_parameter_extremes( aif_0, aif_1 );

		// set up contrast dynamics for healthy liver
		TissueParameter healthy_tissue_0 = mr_cont_gen.get_petmr_tissue_parameter( healthy_tissue_dynamic_labels[0] );
		TissueParameter healthy_tissue_1 = healthy_tissue_0;

		healthy_tissue_0.mr_tissue_.t1_miliseconds_ = t1_healthy_tissue_0;

		healthy_tissue_contrast.set_parameter_extremes( healthy_tissue_0, healthy_tissue_1 );

		// set up contrast dynamics for lesion
		TissueParameter lesion_tissue_0 = mr_cont_gen.get_petmr_tissue_parameter( lesion_dynamic_labels[0] );
		TissueParameter lesion_tissue_1 = lesion_tissue_0;

		lesion_tissue_0.mr_tissue_.t1_miliseconds_ = t1_lesion_0;

		lesion_contrast.set_parameter_extremes( lesion_tissue_0, lesion_tissue_1 );

		// read in the contrast signal 
		SignalContainer aif_dyn_signal = data_io::read_surrogate_signal(filename_contrast_timepoints, filename_aif_t1_);
		SignalContainer healthy_dyn_tissue_signal = data_io::read_surrogate_signal(filename_contrast_timepoints, filename_healthy_tissue_t1);
		SignalContainer lesion_dyn_signal = data_io::read_surrogate_signal(filename_contrast_timepoints, filename_lesion_t1_);


		aif_contrast.set_dyn_signal( aif_dyn_signal );
		healthy_tissue_contrast.set_dyn_signal( healthy_dyn_tissue_signal );
		lesion_contrast.set_dyn_signal( lesion_dyn_signal );


	 	aif_contrast.bin_mr_acquisitions( all_acquis );
		healthy_tissue_contrast.bin_mr_acquisitions( all_acquis );
		lesion_contrast.bin_mr_acquisitions( all_acquis );

		if( num_contrast_states > 0)
		{
			mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (aif_contrast) );
			mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (healthy_tissue_contrast) );
			mr_dyn_sim.add_dynamic( std::make_shared<MRContrastDynamic> (lesion_contrast) );
		}
		
		// ####################################################################################################

		if( simulate_data )
		{
			clock_t t;
			t = clock();
			mr_dyn_sim.simulate_dynamics();
			t = clock() - t;
			
			std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;
			std::stringstream outname_stream;

			outname_stream << "output_grpe_dce_simulation_num_motion_states_" << num_simul_motion_dyn <<"_num_contrast_states_" << num_contrast_states;
			std::string const filename_dce_output = output_path + outname_stream.str() + ".h5";
			mr_dyn_sim.write_simulation_results( filename_dce_output );

		}
		if( store_gt_mvfs && num_simul_motion_dyn > 0 )
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




bool test_pet_dynsim::test_set_template_acquisition_data()
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
		pet_dyn_sim.set_output_filename_prefix("/media/sf_SharedFolder/CCPPETMR/HackathonSimulations/PET/SimulationData");
		
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
		// auto resp_motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		// pet_dyn_sim.add_dynamic( std::make_shared<PETMotionDynamic> (resp_dyn) );
		// resp_dyn.set_displacement_fields( resp_motion_fields, false );
		
		card_dyn.bin_total_time_interval( total_time );
		auto card_motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		card_dyn.set_displacement_fields( card_motion_fields, true );
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

		bool const do_cardiac_sim = false;
		bool const simulate_data = true;
		bool const store_gt_mvfs = false;

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
		
		int const num_sim_motion_states = 4;

		std::cout << "WARNING: NOISE IS STRONGLY SUPPRESSED" << std::endl;
		float const noise_suppression = 1000 * 1000;
		float tot_time_ms =  noise_suppression * 30 * 60 * 1000; // in case there is no motion simulation use this acqu time

		if( num_sim_motion_states > 1)
		{
			PETMotionDynamic motion_dyn( num_sim_motion_states );

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

			SignalContainer motion_signal = data_io::read_surrogate_signal( path_time_pts.str(), path_sig_pts.str());

			// have constant signal
			// std::cout << "WARNING: CONSTANT SIGNAL ASSUMED" << std::endl;
			// for(int i=0; i<motion_signal.size(); i++)
			// 	motion_signal[i].second = 0.0; 
			

			auto first_card_pt = motion_signal[0];
			auto last_card_pt = motion_signal[ motion_signal.size()-1 ];

			float min_time_card_ms = first_card_pt.first;


			for( size_t i=0; i<motion_signal.size(); i++)
			{
				auto curr_sig_pt = motion_signal[i];	
				curr_sig_pt.first = curr_sig_pt.first - min_time_card_ms;
				motion_signal[i] = curr_sig_pt;
			}

		 	motion_dyn.set_dyn_signal( motion_signal );

			tot_time_ms = last_card_pt.first - first_card_pt.first;
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
			clock_t t;
			t = clock();
		
			std::cout << "Simulating Data" << std::endl;
			pet_dyn_sim.simulate_dynamics( tot_time_ms );
			std::cout << "Finished Simulating Data" << std::endl;

			t = clock() - t;
			std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;

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



bool test_pet_dynsim::test_5d_pet_acquisition()
{

	try
	{
		bool const simulate_data = true;
		bool const store_gt_mvfs = false;

		std::string const input_path = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Input/";
		std::string const signal_path = input_path + "/SurrogateSignals/";

		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		PETContrastGenerator pet_cont_gen( segmentation_labels, XML_XCAT_PATH);

		pet_cont_gen.set_template_image_from_file( PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH );

		PETDynamicSimulation pet_dyn_sim( pet_cont_gen );

		// fix name of output
		stringstream prefix_stream;
		prefix_stream << "pet_dyn_5D_"<< "cardiorespiratory";
		
		pet_dyn_sim.set_output_filename_prefix(prefix_stream.str());
		
		pet_dyn_sim.set_filename_rawdata( PET_TEMPLATE_ACQUISITION_DATA_PATH );
		pet_dyn_sim.set_template_image_data( PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH );
		
		int const num_sim_card_states = 8;
		int const num_sim_resp_states = 8;


		float tot_time_ms = 30*60*1000;

		if( num_sim_card_states*num_sim_resp_states > 0 )
		{
			PETMotionDynamic card_dyn( num_sim_card_states ), resp_dyn( num_sim_resp_states );

			std::string fname_timepts, fname_signalpts;

			// read cardiac signal
			fname_timepts = signal_path + "card_time";
			fname_signalpts = signal_path + "card_signal";
			SignalContainer card_signal = data_io::read_surrogate_signal(fname_timepts, fname_signalpts);
			card_dyn.set_dyn_signal( card_signal );
	
			float tot_time_card_ms = aux_test::prep_pet_motion_dyn(card_dyn, card_signal);

			// read cardiac signal
			fname_timepts  = signal_path + "resp_time";
			fname_signalpts = signal_path + "resp_signal";
			SignalContainer resp_signal = data_io::read_surrogate_signal(fname_timepts, fname_signalpts);
			
			float tot_time_resp_ms = aux_test::prep_pet_motion_dyn(resp_dyn, resp_signal);

			std::cout << epiph(tot_time_card_ms) << std::endl;
			std::cout << epiph(tot_time_resp_ms) << std::endl;

			tot_time_ms = 0.5*(tot_time_card_ms+tot_time_resp_ms); // replace total time 


			// read motion fields
			auto card_mvfs = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			card_dyn.set_displacement_fields( card_mvfs, true );
			
			auto resp_mvfs = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			resp_dyn.set_displacement_fields( resp_mvfs, false );
			
			// add the dynamics
			pet_dyn_sim.add_dynamic( std::make_shared<PETMotionDynamic> (card_dyn) );
			pet_dyn_sim.add_dynamic( std::make_shared<PETMotionDynamic> (resp_dyn) );

		}
		if( simulate_data )
		{

			clock_t t;
			t = clock();
		
			std::cout << "Simulating Data" << std::endl;
			pet_dyn_sim.simulate_dynamics( tot_time_ms );
			std::cout << "Finished Simulating Data" << std::endl;

			t = clock() - t;
			std::cout << " TIME FOR SIMULATION: " << (float)t/CLOCKS_PER_SEC/60.f << " MINUTES." <<std::endl;
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
