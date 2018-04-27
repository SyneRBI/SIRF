/* ================================================

Author: Johannes Mayer
Date: 2018.04.03
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include <string>
#include <vector>

#include <ismrmrd/ismrmrd.h>

#include "../auxiliary_input_output.h"

#include "tests_auxiliary_input_output.h"



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
	name_stream << "/media/sf_SharedFiles/test_binary_writer_" << Nx << "x" << Ny << "x" << Nz;

	data_io::write_raw<complex_float_t>(name_stream.str(), dummy_data.begin(), dummy_data.getNumberOfElements());

}