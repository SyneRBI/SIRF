/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "tests_dynamics.h"

#include "auxiliary_input_output.h"
#include "auxiliary_testing_functions.h"
#include "sirf/cReg/NiftiImageData3D.h"

using namespace sirf;

using std::cout;
using std::endl;

bool test_dynamic::test_is_in_bin( void )
{

	bool test_succesful;
	SignalAxisType test_signal = 0.15;

	SignalAxisType bin_min = 0.0;
	SignalAxisType bin_ctr = 0.1;
	SignalAxisType bin_max = 0.2;

	SignalBin bin{bin_min, bin_ctr, bin_max};

	test_succesful = is_in_bin( test_signal , bin );


	std::get<0>(bin) = 0.8;
	std::get<1>(bin) = 0.0;
	std::get<2>(bin) = 0.2;

	test_succesful *= is_in_bin(test_signal, bin);


	std::get<0>(bin) = 0.0;
	std::get<1>(bin) = 0.0;
	std::get<2>(bin) = 0.0;

	test_succesful *= !(is_in_bin(test_signal, bin));

	return test_succesful;

}


bool test_dynamic::test_intersect_mr_acquisition_data( void )
{

	AcquisitionsVector one_vec, other_vec;
	auto hdr = aux_test::get_serialized_mock_ismrmrd_header();
	
	one_vec.set_acquisitions_info( hdr );
	other_vec.set_acquisitions_info( hdr );
 

	ISMRMRD::Acquisition acq;
	ISMRMRD::AcquisitionHeader acq_hdr = acq.getHead();

	uint32_t const one_start_counter = 1;
	uint32_t const one_end_counter = 10;

	uint32_t const other_start_counter = 6;
	uint32_t const other_end_counter = 13;

	
	for( uint32_t i=one_start_counter; i<=one_end_counter; i++)
	{	
		acq_hdr.scan_counter = i;
		acq.setHead( acq_hdr );
		one_vec.append_acquisition(acq);
	}

	for( uint32_t i=other_start_counter; i<=other_end_counter; i++)
	{	
		acq_hdr.scan_counter = i;
		acq.setHead( acq_hdr );
		other_vec.append_acquisition(acq);
	}
	
	AcquisitionsVector intersec_vec = intersect_mr_acquisition_data(one_vec, other_vec);

	int32_t num_overlap = one_end_counter - other_start_counter;
	num_overlap = (num_overlap>=0) ? num_overlap+1 : 0;


	bool test_succesful = (num_overlap == intersec_vec.items());
	
	cout << epiph(num_overlap) << endl;
	cout << epiph(intersec_vec.items()) << endl;


	return test_succesful;
}


bool test_dynamic::test_linear_interpolate_signal( )
{

	try
	{
		SignalContainer mock_signal = aux_test::get_mock_motion_signal();

		aDynamic dyn;
		dyn.set_dyn_signal(mock_signal);

		
		size_t num_repetitions_for_speed_estimate = 256*256*3;

		for( size_t rep=0; rep<num_repetitions_for_speed_estimate; rep++)
		{
			TimeAxisType time_point = 0.55;
			SignalAxisType interpol_signal = dyn.linear_interpolate_signal( time_point );
		}
		TimeAxisType time_point = 0.55;
		SignalAxisType interpol_signal = dyn.linear_interpolate_signal( time_point );

		cout << epiph ( interpol_signal ) <<endl; 

		time_point = -1;
		interpol_signal = dyn.linear_interpolate_signal( time_point );
		cout << epiph ( interpol_signal ) <<endl;

		time_point = 15;
		interpol_signal = dyn.linear_interpolate_signal( time_point );
		cout << epiph ( interpol_signal ) <<endl;	

		return true;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}
}


bool test_dynamic::test_get_set_bins()
{
	try
	{
		bool test_succesful = true;

		int const num_bins = 10;
		aDynamic dyn(num_bins);

		std::vector< SignalBin > all_bins = dyn.get_bins();

		std::cout << "#### REGULAR DYNAMIC ####" << std::endl;
		for( int i=0; i<all_bins.size(); i++ )
		{
			cout << epiph(std::get<0> (all_bins[i] )) << endl;
			cout << epiph(std::get<1> (all_bins[i] )) << endl;
			cout << epiph(std::get<2> (all_bins[i] )) << endl;
		}

		test_succesful = ( all_bins.size() == num_bins );


		ContrastDynamic cont_dyn(num_bins);

		all_bins = cont_dyn.get_bins();

		std::cout << "#### CONTRAST DYNAMIC ####" << std::endl;
		for( int i=0; i<all_bins.size(); i++ )
		{
			cout << epiph(std::get<0> (all_bins[i] )) << endl;
			cout << epiph(std::get<1> (all_bins[i] )) << endl;
			cout << epiph(std::get<2> (all_bins[i] )) << endl;
		}

		test_succesful = ( all_bins.size() == num_bins+1 );
		
		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}
}



bool test_dynamic::test_bin_mr_acquisitions()
{
	try
	{
		bool test_succesful = true;

		sirf::AcquisitionsVector acq_vec;
		acq_vec.read( std::string(ISMRMRD_H5_TEST_PATH ));

		int const num_bins = 10;
		MRMotionDynamic motion_dyn(num_bins);

		SignalContainer mock_signal = aux_test::get_generic_cardiac_signal(acq_vec);
		motion_dyn.set_dyn_signal( mock_signal );

		motion_dyn.bin_mr_acquisitions( acq_vec );
		auto motion_binned_acquis = motion_dyn.get_binned_mr_acquisitions();

		test_succesful *= (motion_binned_acquis.size() == num_bins);

		std::cout << "######## Motion Dynamic: ########"<< std::endl;
		size_t num_tot_acquis = 0;
		for( int i=0; i<num_bins; i++)
		{
			num_tot_acquis += motion_binned_acquis[i].items();
			cout << epiph( motion_binned_acquis[i].items() ) << endl;
		}

		cout << epiph(num_tot_acquis) <<endl;
		cout << epiph(acq_vec.items()) <<endl;

		test_succesful *= ( num_tot_acquis == acq_vec.items() );


		
		std::cout << "######## Contrast Dynamic: ########"<< std::endl;

		
		MRContrastDynamic mr_cont_dyn( num_bins );
		std::cout << epiph( mr_cont_dyn.get_num_simul_states() ) << std::endl;

		mr_cont_dyn.set_dyn_signal( mock_signal );
		mr_cont_dyn.bin_mr_acquisitions( acq_vec );

		auto contrast_binned_acquis = mr_cont_dyn.get_binned_mr_acquisitions();
		test_succesful *= (contrast_binned_acquis.size() == num_bins+1);

		std::cout << epiph( contrast_binned_acquis.size() ) << std::endl;

		num_tot_acquis = 0;
		for( int i=0; i<contrast_binned_acquis.size(); i++)
		{
			num_tot_acquis += contrast_binned_acquis[i].items();
			cout << epiph( contrast_binned_acquis[i].items() ) << endl;
		}

		cout << epiph(num_tot_acquis) <<endl;
		cout << epiph(acq_vec.items()) <<endl;

		test_succesful *= ( num_tot_acquis == acq_vec.items() );
		
		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}


}


bool test_dynamic::test_motion_dynamic_counter()
{
	try
	{
		bool test_succesful = true;

		MotionDynamic first_dyn, second_dyn, third_dyn;

		
		int const first_counter = first_dyn.get_which_motion_dynamic_am_i();
		int const second_counter = second_dyn.get_which_motion_dynamic_am_i();
		int const third_counter = third_dyn.get_which_motion_dynamic_am_i();

		int const total_counter = first_dyn.get_num_total_motion_dynamics();

		cout << epiph(first_counter) << endl;
		cout << epiph(second_counter) << endl;
		cout << epiph(third_counter) << endl;
		cout << epiph(total_counter) << endl;

		test_succesful *= (first_counter == 0);
		test_succesful *= (second_counter == 1);
		test_succesful *= (third_counter == 2);
		test_succesful *= (total_counter == 3);

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}
}




bool test_dynamic::test_motion_dynamic_temp_folder_setup( )
{

	try
	{
		bool test_succesful = true;

		MotionDynamic first_dyn;
		cout << epiph( first_dyn.get_temp_folder_name() )<< endl;
		// test_succesful *= first_dyn.make_temp_folder();
		// test_succesful *= first_dyn.delete_temp_folder();

		// test_succesful *= !(first_dyn.delete_temp_folder());

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}

}

bool test_dynamic::test_motion_dynamic_set_motion_fields()
{
try
	{
		bool test_succesful = true;

		auto resp_mvfs = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );

		int const num_simul_bins = 12;
		MotionDynamic motion_dyn(num_simul_bins);
		motion_dyn.set_displacement_fields(resp_mvfs);

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}

}

bool test_dynamic::test_mvf_vs_pet_img_quarternions( void )
{
	try
	{
		bool test_succesful = true;

		auto pet_img = STIRImageData(PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH);
		NiftiImageData3D<float> pet_sirf_img( pet_img );
		auto pet_img_data = pet_sirf_img.get_raw_nifti_sptr();

		float const img_off_x = pet_img_data->qoffset_x;
		float const img_off_y = pet_img_data->qoffset_y;
		float const img_off_z = pet_img_data->qoffset_z;
		 

		float const img_quart_b = pet_img_data->quatern_b;
		float const img_quart_c = pet_img_data->quatern_c;
		float const img_quart_d = pet_img_data->quatern_d;
		float const img_quart_ac = pet_img_data->qfac;

		float const img_slope = pet_img_data->scl_slope;
		float const img_inter = pet_img_data->scl_inter;

		float const img_dx = pet_img_data->dx;
	    float const img_dy = pet_img_data->dy;
	    float const img_dz = pet_img_data->dz;
	    float const img_dt = pet_img_data->dt;
	    float const img_du = pet_img_data->du;
	    float const img_dv = pet_img_data->dv;
	    float const img_dw = pet_img_data->dw;

		cout << epiph( img_off_x )<< endl;
		cout << epiph( img_off_y )<< endl;
		cout << epiph( img_off_z )<< endl;

		cout << epiph( img_quart_b ) << endl;
		cout << epiph( img_quart_c ) << endl;
		cout << epiph( img_quart_d ) << endl;
		cout << epiph( img_quart_ac) << endl;
        cout << epiph( img_slope ) << endl;
		cout << epiph( img_inter ) << endl;


		cout << epiph( img_dx ) << endl;
		cout << epiph( img_dy ) << endl;
		cout << epiph( img_dz ) << endl;
		cout << epiph( img_dt ) << endl;
		cout << epiph( img_du ) << endl;
		cout << epiph( img_dv ) << endl;
		cout << epiph( img_dw ) << endl;

		auto resp_mvfs = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );

		int const num_simul_bins = 1;
		PETMotionDynamic motion_dyn(num_simul_bins);
		motion_dyn.set_displacement_fields(resp_mvfs);
		motion_dyn.prep_displacements_fields();
		motion_dyn.align_motion_fields_with_image(pet_img);


		auto some_mvf = motion_dyn.get_interpolated_displacement_field( 0.f );
		auto mvf_img_data = some_mvf.get_raw_nifti_sptr();

		float const mvf_off_x = mvf_img_data->qoffset_x;
		float const mvf_off_y = mvf_img_data->qoffset_y;
		float const mvf_off_z = mvf_img_data->qoffset_z;

		float const mvf_quart_b = mvf_img_data->quatern_b;
		float const mvf_quart_c = mvf_img_data->quatern_c;
		float const mvf_quart_d = mvf_img_data->quatern_d;
		float const mvf_quart_ac = mvf_img_data->qfac;

		float const mvf_slope = mvf_img_data->scl_slope;
		float const mvf_inter = mvf_img_data->scl_inter;

		float const mvf_dx = mvf_img_data->dx;
	    float const mvf_dy = mvf_img_data->dy;
	    float const mvf_dz = mvf_img_data->dz;
	    float const mvf_dt = mvf_img_data->dt;
	    float const mvf_du = mvf_img_data->du;
	    float const mvf_dv = mvf_img_data->dv;
	    float const mvf_dw = mvf_img_data->dw;

		cout << epiph( mvf_off_x)<< endl;
		cout << epiph( mvf_off_y)<< endl;
		cout << epiph( mvf_off_z)<< endl;

		cout << epiph(mvf_quart_b ) << endl;
		cout << epiph(mvf_quart_c ) << endl;
		cout << epiph(mvf_quart_d ) << endl;
		cout << epiph(mvf_quart_ac) << endl;

		cout << epiph( mvf_slope ) << endl;
		cout << epiph( mvf_inter ) << endl;
		
		cout << epiph( mvf_dx ) << endl;
		cout << epiph( mvf_dy ) << endl;
		cout << epiph( mvf_dz ) << endl;
		cout << epiph( mvf_dt ) << endl;
		cout << epiph( mvf_du ) << endl;
		cout << epiph( mvf_dv ) << endl;
		cout << epiph( mvf_dw ) << endl;

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}


}

bool test_dynamic::test_motion_dynamic_prep_motion_fields()
{
try
	{
		bool test_succesful = true;

		auto resp_mvfs = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );

		int const num_simul_bins = 12;
		MotionDynamic motion_dyn(num_simul_bins);
		motion_dyn.set_displacement_fields(resp_mvfs);
		motion_dyn.prep_displacements_fields();

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}

}



bool test_dynamic::test_motion_dynamic_temp_interpolate_dvfs( void )
{
	try
	{
		bool test_succesful = true;

		auto resp_mvfs = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );

		int const num_simul_bins = 12;
		MotionDynamic motion_dyn(num_simul_bins);
		motion_dyn.set_displacement_fields(resp_mvfs);
		motion_dyn.prep_displacements_fields();

		SignalAxisType motion_signal = 0.6;

		auto interpolated_dvf = motion_dyn.get_interpolated_displacement_field(motion_signal);

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}
}



bool test_dynamic::test_mr_contrast_motion_dyn_get_num_simul_states( void )
{
	try
	{
		bool test_succesful = true;

		int const num_con_states = 3;
		int const num_motion_states = 4;

		MRContrastDynamic cont_dyn( num_con_states );
		MRMotionDynamic motion_dyn( num_motion_states );

		cout << epiph (cont_dyn.get_num_simul_states()) << endl;
		cout << epiph (motion_dyn.get_num_simul_states()) << endl;

		test_succesful *= (cont_dyn.get_num_simul_states() == num_con_states);
		test_succesful *= (motion_dyn.get_num_simul_states() == num_motion_states);

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}
}



bool test_dynamic::test_bin_pet_time_interval( void )
{
	try
	{
		bool test_succesful = true;
		// int const num_sig_pts = 4;
		// for(int i=0; i<num_sig_pts; i++)
		// {
		// 	SignalPoint sp;
		// 	sp.first = i*100 + i*i*sqrt(2);
		// 	sp.second = (SignalAxisType)(i%2);
		// 	signal_cont.push_back( sp );
		// }

		SignalContainer card_sig = data_io::read_surrogate_signal( std::string(TIME_POINTS_CARDIAC_PATH), std::string(CARDIAC_SIGNAL_PATH));
		
		int const num_simul_bins = 10;
		aPETDynamic pet_dyn(num_simul_bins);

	
		auto first_pt = card_sig[0];
		auto last_pt = card_sig[ card_sig.size()-1 ];

		float min_time_ms = first_pt.first;
		float tot_time_ms = last_pt.first - first_pt.first;

		std::cout << "total time ms: " << tot_time_ms <<std::endl;

		for( size_t i=0; i<card_sig.size(); i++)
		{
			auto curr_sig_pt = card_sig[i];	
			curr_sig_pt.first = 25000 * (curr_sig_pt.first - min_time_ms)/tot_time_ms;
			std::cout << curr_sig_pt.first <<std::endl;
			card_sig[i] = curr_sig_pt;
		}

	 	first_pt = card_sig[0];
		last_pt = card_sig[ card_sig.size()-1 ];

		min_time_ms = first_pt.first;
		tot_time_ms = last_pt.first - first_pt.first;

		std::cout << "total time ms: " << tot_time_ms <<std::endl;
	 	pet_dyn.set_dyn_signal( card_sig );

	 	pet_dyn.bin_total_time_interval( TimeBin(0,tot_time_ms) );
		TimeBin total_time(0, 25);

		pet_dyn.bin_total_time_interval( total_time );

		for(int i=0; i<num_simul_bins; i++)
			cout <<"Time spent in bin " <<i << " = " << pet_dyn.get_time_spent_in_bin(i) <<endl;




		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}
}