/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "tests_dynamics.h"

#include "auxiliary_input_output.h"
#include "auxiliary_testing_functions.h"

using namespace sirf;

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
	
	std::cout << epiph(num_overlap) << std::endl;
	std::cout << epiph(intersec_vec.items()) << std::endl;


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

		std::cout << epiph ( interpol_signal ) <<std::endl; 

		time_point = -1;
		interpol_signal = dyn.linear_interpolate_signal( time_point );
		std::cout << epiph ( interpol_signal ) <<std::endl;

		time_point = 15;
		interpol_signal = dyn.linear_interpolate_signal( time_point );
		std::cout << epiph ( interpol_signal ) <<std::endl;	

		return true;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
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


		for( int i=0; i<all_bins.size(); i++ )
		{
			std::cout << epiph(std::get<0> (all_bins[i] )) << std::endl;
			std::cout << epiph(std::get<1> (all_bins[i] )) << std::endl;
			std::cout << epiph(std::get<2> (all_bins[i] )) << std::endl;
		}

		test_succesful = ( all_bins.size() == num_bins );
		
		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}



bool test_dynamic::test_bin_mr_acquisitions()
{
	try
	{
		bool test_succesful = true;

		auto acq_vec = mr_io::read_ismrmrd_acquisitions( ISMRMRD_H5_TEST_PATH );

		int const num_bins = 10;
		aDynamic dyn(num_bins);

		SignalContainer mock_signal = aux_test::get_mock_motion_signal(acq_vec);
		dyn.set_dyn_signal( mock_signal );

		dyn.bin_mr_acquisitions( acq_vec );
		auto binned_acquis = dyn.get_binned_mr_acquisitions();

		test_succesful *= (binned_acquis.size() == num_bins);

		size_t num_tot_acquis = 0;
		for( int i=0; i<num_bins; i++)
		{
			num_tot_acquis += binned_acquis[i].items();
			std::cout << epiph( binned_acquis[i].items() ) << std::endl;
		}

		std::cout << epiph(num_tot_acquis) <<std::endl;
		std::cout << epiph(acq_vec.items()) <<std::endl;


		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
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

		std::cout << epiph(first_counter) << std::endl;
		std::cout << epiph(second_counter) << std::endl;
		std::cout << epiph(third_counter) << std::endl;
		std::cout << epiph(total_counter) << std::endl;

		test_succesful *= (first_counter == 0);
		test_succesful *= (second_counter == 1);
		test_succesful *= (third_counter == 2);
		test_succesful *= (total_counter == 3);

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}




