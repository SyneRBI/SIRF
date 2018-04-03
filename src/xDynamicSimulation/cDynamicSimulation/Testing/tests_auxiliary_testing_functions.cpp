/* ================================================

Author: Johannes Mayer
Date: 2018.04.03
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tests_auxiliary_testing_functions.h"


void test_aux_test_funs::test_write_ndarray_to_raw( void )
{

	size_t Nx = 3;
	size_t Ny = 3; 
	size_t Nz = 3;
	size_t Ne = 4;

	std::vector< size_t > data_size = {Nx, Ny, Nz, Ne};


	ISMRMRD::NDArray< complex_float_t > dummy_data;
	dummy_data.resize(data_size);

	for( int nz=0; nz<Nz; nz++)
	for( int ny=0; ny<Ny; ny++)
	for( int nx=0; nx<Nx; nx++)
	for( int ne=0; ne<Ne; ne++)
	{
		dummy_data(nx,ny,nz,ne) = std::complex<float>(ne, 0);
	}


	std::string output_name = "/media/sf_SharedFiles/test_binary_writer";

	aux_test::write_ndarray_to_binary(output_name, dummy_data);

}