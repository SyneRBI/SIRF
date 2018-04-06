/*  ###############################################
author: johannes mayer
date: 15. March 2018


##################*/
#pragma once


#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>



class aFullySampledFFT{

public:

	aFullySampledFFT(ISMRMRD::IsmrmrdHeader hdr)
	{
		this->hdr_ = hdr;
	}

	virtual void SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data) = 0;

protected:

	ISMRMRD::NDArray<complex_float_t> k_data_;
	ISMRMRD::IsmrmrdHeader hdr_;

};


class FullySampledCartesianFFT: public aFullySampledFFT{

public:
	FullySampledCartesianFFT(ISMRMRD::IsmrmrdHeader hdr);

};
