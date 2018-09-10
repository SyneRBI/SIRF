/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <ismrmrd/ismrmrd.h>

#include <gadgetron/hoNDArray.h>
#include <gadgetron/ho2DArray.h>
#include <gadgetron/vector_td.h>


// #include <ismrmrd/xml.h>

typedef float TrajectoryPrecision;
typedef ISMRMRD::NDArray< TrajectoryPrecision > TrajectoryContainer;

typedef Gadgetron::floatd2 TrajectoryType2D;





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

template < typename TrajType >
class aTrajectoryPreparation 
{

public:

	std::vector< size_t > get_traj_dims()
	{
		return traj_dims_;
	}

	Gadgetron::hoNDArray< TrajType > get_formatted_trajectory( void ) 
	{
		return this->traj_;
	};

	template< typename OutputDataType >
	Gadgetron::ho2DArray< OutputDataType > get_formatted_output_container( void )
	{
		if(this->traj_dims_.size() == 0)
			throw std::runtime_error("The trajectory has not been set. Please do so before calling this function.");
		Gadgetron::ho2DArray< OutputDataType > output_container( &(this->traj_dims_) );

		output_container.fill( OutputDataType(0) );


		return output_container;
	}

	template< typename OutputDataType >
	Gadgetron::hoNDArray< OutputDataType > get_formatted_identity_dcf( void )
	{
		if(this->traj_dims_.size() == 0)
			throw std::runtime_error("The trajectory has not been set. Please do so before calling this function.");
		Gadgetron::hoNDArray< OutputDataType > identity_dcf( this->traj_dims_ );

		identity_dcf.fill( OutputDataType(1.f) );

		return identity_dcf;
	}


protected:

	std::vector< size_t > traj_dims_;
	Gadgetron::hoNDArray< TrajType > traj_;	


};


class RPETrajectoryPreparation: public aTrajectoryPreparation< TrajectoryType2D >{

public:
	void set_and_check_trajectory( ISMRMRD::NDArray< TrajectoryPrecision > trajectory);

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

	RPETrajectoryPreparation traj_prep_;
	
};




