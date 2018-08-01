/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "tests_dynamics.h"

#include "auxiliary_testing_functions.h"




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












