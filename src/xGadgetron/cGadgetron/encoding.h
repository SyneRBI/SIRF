/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

class aFullySampledFFT{

public:

	aFullySampledFFT()
	{}

	virtual ISMRMRD::NDArray<complex_float_t> get_k_data( void );

	virtual void SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data ) = 0;

protected:

	ISMRMRD::NDArray<complex_float_t> k_data_;
	
};


class FullySampledCartesianFFT: public aFullySampledFFT{

public:
	FullySampledCartesianFFT();
	void SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data);

};
