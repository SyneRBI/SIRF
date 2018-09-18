/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "tests_encoding.h"

#include <string>

#include <ismrmrd/ismrmrd.h>

#include <gadgetron/vector_td.h>

#include "gadgetron_data_containers.h"

#include "auxiliary_input_output.h"
#include "auxiliary_testing_functions.h"


using ISMRMRD::NDArray;


bool test_enc::test_cube_input()
{

	NDArray<complex_float_t> arr = aux_test::get_mock_ndarray_with_cube();

	std::string output_name = std::string( SHARED_FOLDER_PATH ) + "test_enc_mock_array" ;

	size_t num_elements = arr.getNumberOfElements();

	std::vector<float> arr_abs, arr_arg;
	arr_abs.resize(num_elements);
	arr_arg.resize(num_elements);

	for( size_t i=0; i<num_elements; i++)
	{
		arr_abs[i] = std::abs( *(arr.begin() + i));
		arr_arg[i] = std::arg( *(arr.begin() + i));
	}

	data_io::write_raw<float>(output_name + "_abs_64x64x64", &arr_abs[0], arr_abs.size());
	data_io::write_raw<float>(output_name + "_arg_64x64x64", &arr_arg[0], arr_arg.size());

	return true;
}

bool CartesianEncodingTester::test_sample_fourier_space()
{
	FullySampledCartesianFFT cart_fft;

	NDArray<complex_float_t> i_dat = aux_test::get_mock_ndarray_with_cube();
	
	cart_fft.SampleFourierSpace( i_dat);

	NDArray<complex_float_t> k_dat = cart_fft.get_k_data();

	size_t num_elements = k_dat.getNumberOfElements();

	std::vector<float> k_dat_abs;

	k_dat_abs.resize(num_elements);

	for( size_t i=0; i<num_elements; i++)
		k_dat_abs[i] = std::abs( *(k_dat.begin() + i) );

	
	std::string output_name =  std::string(SHARED_FOLDER_PATH)+ "test_cart_enc_k_data" ;
	data_io::write_raw<float>(output_name + "_abs_64x64x64", &k_dat_abs[0], k_dat_abs.size());

	return true;
}



bool RPETester::test_sample_fourier_space( void )
{
	try
	{

		
		NDArray<complex_float_t> i_dat = aux_test::get_mock_ndarray_with_cube();
		
		auto img_dims = i_dat.getDims();

		size_t const NRad = 64;
		size_t const NAng = 64;
		
		sirf::RPETrajectoryContainer rpe_traj = aux_test::get_mock_radial_trajectory(NRad, NAng);

		auto radial_traj = rpe_traj.get_trajectory();

		RadialPhaseEncodingFFT rpe_fft;
		
		rpe_fft.set_trajectory( radial_traj );
		
		rpe_fft.SampleFourierSpace( i_dat);
		
		NDArray<complex_float_t> k_dat = rpe_fft.get_k_data();

		size_t num_elements = k_dat.getNumberOfElements();

		std::vector<float> k_dat_abs;

		k_dat_abs.resize(num_elements);

		for( size_t i=0; i<num_elements; i++)
			k_dat_abs[i] = std::abs( *(k_dat.begin() + i) );

		std::stringstream output_name;
		output_name << std::string(SHARED_FOLDER_PATH)+ "test_rpe_enc_k_data_" ;
		output_name << img_dims[0] << "x" << NRad << "x"<< NAng;

		data_io::write_raw<float>(output_name.str(), &k_dat_abs[0], k_dat_abs.size());

		return true;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}



bool RPETrajectoryPreparationTester::test_get_set_trajectory()
{
try
	{
		size_t const NRad = 64;
		size_t const NAng = 64;

		
		sirf::RPETrajectoryContainer rpe_traj = aux_test::get_mock_radial_trajectory(NRad, NAng);
		auto radial_traj = rpe_traj.get_trajectory();

		RPETrajectoryPreparation traj_prep;
		traj_prep.set_and_check_trajectory( radial_traj);

		auto formatted_traj = traj_prep.get_formatted_trajectory();

		return true;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}	


bool RPETrajectoryPreparationTester::test_get_result_container()
{
try
	{
		size_t const NRad = 64;
		size_t const NAng = 64;

		sirf::RPETrajectoryContainer rpe_traj = aux_test::get_mock_radial_trajectory(NRad, NAng);
		auto radial_traj = rpe_traj.get_trajectory();

		RPETrajectoryPreparation traj_prep;
		traj_prep.set_and_check_trajectory( radial_traj);

		auto formatted_output = traj_prep.get_formatted_output_container< complex_float_t >();

		return true;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}	