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
