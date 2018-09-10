/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <ismrmrd/ismrmrd.h>
// #include <ismrmrd/xml.h>

// typedef std::vector<float> TrajectoryType;
typedef float TrajectoryType;
typedef ISMRMRD::NDArray< TrajectoryType > TrajectoryContainer;



template <typename T>
std::vector<size_t> data_dims_from_ndarray(ISMRMRD::NDArray< T > data)
{
	std::vector<size_t> data_dims( data.getDims(), data.getDims()+ISMRMRD::ISMRMRD_NDARRAY_MAXDIM );
	for( int i=0; i<ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i++)
	{
		if( data_dims[i] == 0)
			data_dims[i] = 1;
	}
	return data_dims;
};




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

	void set_trajectory(TrajectoryContainer &traj);
	void SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data );

private:

	void prep_trajctory( void );

	TrajectoryContainer traj_;

};
