/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "tests_dynamics.h"

#include "auxiliary_input_output.h"
#include "auxiliary_testing_functions.h"


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

		test_succesful = (binned_acquis.size() == num_bins);

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}


}




