/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "encoding.h"

#include <gadgetron/hoNDArray.h>
#include <gadgetron/hoNDFFT.h>

using Gadgetron::hoNDFFT;

ISMRMRD::NDArray<complex_float_t> aFullySampledFFT::get_k_data( void )
{
	return this->k_data_;	
}


FullySampledCartesianFFT::FullySampledCartesianFFT():
aFullySampledFFT()
{
}

void FullySampledCartesianFFT::SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data)
{

	size_t const num_elements = i_data.getNumberOfElements();
	std::vector<size_t> data_dims( i_data.getDims(), i_data.getDims()+ISMRMRD::ISMRMRD_NDARRAY_MAXDIM );
	for( int i=0; i<ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i++)
		{
			if( data_dims[i] == 0)
				data_dims[i] = 1;
		}

	Gadgetron::hoNDArray< complex_float_t > data_to_be_fftd(data_dims);
	
	for( size_t i=0; i<num_elements; i++)
		*(data_to_be_fftd.begin() + i) = *(i_data.begin() + i);

	// transform first three dims
	hoNDFFT< float >::instance()->fft3c( data_to_be_fftd);
	
	this->k_data_ = ISMRMRD::NDArray<complex_float_t>(data_dims);

	for( size_t i=0; i<num_elements; i++)
	{
		*(this->k_data_.begin() + i) = *(data_to_be_fftd.begin() + i);
	}

}