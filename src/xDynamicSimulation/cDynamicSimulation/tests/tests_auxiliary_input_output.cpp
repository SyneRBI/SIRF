/* ================================================

Author: Johannes Mayer
Date: 2018.04.03
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include <string>
#include <vector>

#include <ismrmrd/ismrmrd.h>

#include "sirf/cDynamicSimulation/auxiliary_input_output.h"
#include "auxiliary_testing_functions.h"

#include "tests_auxiliary_input_output.h"

using namespace sirf;

void test_aux_io::test_write_ndarray_to_raw( void )
{

	size_t Nx = 192;
	size_t Ny = 192; 
	size_t Nz = 192;
	size_t Ne = 3;

	std::vector< size_t > data_size = {Nx, Ny, Nz, Ne};


	ISMRMRD::NDArray< complex_float_t > dummy_data;
	dummy_data.resize(data_size);

	for( int nz=0; nz<Nz; nz++)
	for( int ny=0; ny<Ny; ny++)
	for( int nx=0; nx<Nx; nx++)
	for( int ne=0; ne<Ne; ne++)
	{
		dummy_data(nx,ny,nz,ne) = std::complex<float>(nx*ne, nx*ne);
	}


	std::stringstream name_stream;
	name_stream << SHARED_FOLDER_PATH << "test_binary_writer_" << Nx << "x" << Ny << "x" << Nz;

	data_io::write_raw<complex_float_t>(name_stream.str(), dummy_data.begin(), dummy_data.getNumberOfElements());

}



bool test_aux_io::test_read_acquisitions_vector_number_consistency( void )
{

	size_t const expected_num_acquisitions = 128*128;


	AcquisitionsVector acqu_vec = mr_io::read_ismrmrd_acquisitions(ISMRMRD_H5_TEST_PATH);


	size_t const read_num_acquisitions = acqu_vec.items();	

	std::cout << epiph(expected_num_acquisitions) << std::endl;
	std::cout << epiph(read_num_acquisitions) << std::endl;


	return read_num_acquisitions == expected_num_acquisitions;

}




void test_aux_io::test_write_ismrmrd_image_to_analyze( void )
{
	typedef complex_float_t input_type_mock_object;
	auto img =  aux_test::get_mock_ismrmrd_image_with_cube(  );

	data_io::write_ISMRMRD_Image_to_nii< input_type_mock_object > (ANALYZE_OUTPUT_TESTPATH, img);

}



bool test_aux_io::test_read_single_column_txt( void )
{
	try
	{	
		bool test_successful = true;
		
		// std::string const filename_input = std::string(SHARED_FOLDER_PATH) + "testdata_inputoutput/testfile";
		std::string const filename_input = std::string(SHARED_FOLDER_PATH) + "/PublicationData/Input/SurrogateSignals/card_time";

		std::vector< float > input = data_io::read_single_column_txt<float>(filename_input);
		std::cout << epiph( input.size() ) << std::endl;

		for( size_t i=0; i<input.size(); i++)
		{
			std::cout << epiph( input[i] ) <<std::endl;
		}

		return  test_successful;
	}
	catch(...)
	{
		std::cout << "An unknown exception was caught" << std::endl;
		return false;
	}


}





