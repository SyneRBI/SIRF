/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "tests_encoding.h"

#include <string>

#include <ismrmrd/ismrmrd.h>

#include "auxiliary_input_output.h"
#include "auxiliary_testing_functions.h"


using ISMRMRD::NDArray;


bool test_enc::test_cube_input()
{

	NDArray<complex_float_t> arr = aux_test::get_mock_ndarray_with_cube();

	std::string output_name = SHARED_FOLDER_PATH + "test_enc_mock_array";

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

bool test_cart_enc::test_sample_fourier_space()
{
	ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
	FullySampledCartesianFFT cart_fft;


	NDArray<complex_float_t> i_dat = aux_test::get_mock_ndarray_with_cube();
	
	cart_fft.SampleFourierSpace( i_dat);

	NDArray<complex_float_t> k_dat = cart_fft.get_k_data();

	size_t num_elements = k_dat.getNumberOfElements();

	std::vector<float> k_dat_abs;

	k_dat_abs.resize(num_elements);

	for( size_t i=0; i<num_elements; i++)
		k_dat_abs[i] = std::abs( *(k_dat.begin() + i) );

	
	std::string output_name = SHARED_FOLDER_PATH + "test_cart_enc_k_data";
	data_io::write_raw<float>(output_name + "_abs_64x64x64", &k_dat_abs[0], k_dat_abs.size());

	return true;
}



