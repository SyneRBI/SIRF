/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <ismrmrd/ismrmrd.h>
// #include <ismrmrd/xml.h>

typedef ISMRMRD::NDArray< std::vector<float> > TrajectoryType;



class aCartesianReadoutFFT{

public:

	aCartesianReadoutFFT()
	{}

	virtual ISMRMRD::NDArray<complex_float_t> get_k_data( void );

	virtual void SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data ) = 0;

protected:

	ISMRMRD::NDArray<complex_float_t> k_data_;
	
};


class FullySampledCartesianFFT: public aCartesianReadoutFFT{

public:
	FullySampledCartesianFFT();
	void SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data);

};


class RadialPhaseEncodingFFT: public aCartesianReadoutFFT{

public:

	RadialPhaseEncodingFFT() {};

	void set_trajectory(TrajectoryType &traj);
	void SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data );

private:

	TrajectoryType traj_;

};