/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <utility>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include <gadgetron/hoNDArray.h>
#include <gadgetron/ho2DArray.h>
#include <gadgetron/vector_td.h>

#include "gadgetron_data_containers.h"



namespace sirf
{


typedef float TrajPrecision;
typedef ISMRMRD::NDArray<TrajPrecision> TrajVessel;

typedef Gadgetron::floatd2 TrajectoryType2D;
typedef ISMRMRD::NDArray<complex_float_t> MREncodingDataType;


template <typename T>
std::vector<size_t> data_dims_from_ndarray(ISMRMRD::NDArray< T > &data)
{
	std::vector<size_t> data_dims( data.getDims(), data.getDims()+ISMRMRD::ISMRMRD_NDARRAY_MAXDIM );
	for( int i=0; i<ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i++)
	{
		if( data_dims[i] == 0)
			data_dims[i] = 1;
	}
	return data_dims;
};


/*
\ingroup Gadgetron Data Containers
\brief A trajectory container to pass non-cartesian sampling patterns to an MRAcquisition model.

Acquisition models are stored in ISMRMRD::NDArray<float> containers.
Additional functionality is provided to overwrite ISMRMRDHeader information.

*/

class aTrajectoryContainer{

public:

	aTrajectoryContainer()
	{
		this->traj_.resize( std::vector<size_t>{0} );	
	}
	
	void set_header(ISMRMRD::IsmrmrdHeader hdr);
	void set_trajectory( TrajVessel trajectory );
	
	std::string get_traj_type( void );
	TrajVessel get_trajectory( void );

	void overwrite_ismrmrd_trajectory_info(ISMRMRD::IsmrmrdHeader& hdr);
	void overwrite_ismrmrd_trajectory_info(std::string& serialized_header);

	virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& aqu)=0;
	virtual void compute_trajectory()=0;
protected:

	virtual TrajPrecision get_traj_max_abs( void )=0;
	void norm_trajectory( void );

	ISMRMRD::IsmrmrdHeader hdr_;

	TrajVessel traj_;
	std::string traj_type_;

};

class CartesianTrajectoryContainer : public aTrajectoryContainer{

public:
	CartesianTrajectoryContainer() :aTrajectoryContainer()
	{
		traj_type_ = "Cartesian"; 
	}

	void set_acquisition_trajectory(ISMRMRD::Acquisition& aqu){};
	void compute_trajectory() {};
protected:
	virtual TrajPrecision get_traj_max_abs( void ){ return 0;}
};


class RPETrajectoryContainer : public aTrajectoryContainer{

public:
	RPETrajectoryContainer():aTrajectoryContainer()
	{	
		traj_type_ = "RPE"; 
	}
	
	void set_acquisition_trajectory(ISMRMRD::Acquisition& aqu);
	virtual void compute_trajectory( void );

protected:
	TrajPrecision get_traj_max_abs( void );	

};

class RPEInterleavedTrajectoryContainer: public RPETrajectoryContainer{

public:
	void compute_trajectory( void );
};

class RPEInterleavedGoldenCutTrajectoryContainer: public RPETrajectoryContainer{
	
public:
	void compute_trajectory( void );
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
	Gadgetron::hoNDArray< OutputDataType > get_formatted_output_container( void )
	{
		if(this->traj_dims_.size() == 0)
			throw std::runtime_error("The trajectory has not been set. Please do so before calling this function.");
		Gadgetron::hoNDArray< OutputDataType > output_container( &(this->traj_dims_) );

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
	void set_and_check_trajectory( TrajVessel& trajectory );
	void set_kspace_subset( sirf::MRAcquisitionData &ad );

	void get_traj_index_pair( size_t const idx ) const;
	
protected:
	std::vector< std::pair< size_t, size_t > > lut_idx_to_traj_;

};


class aCartesianReadoutFFT{

public:

	aCartesianReadoutFFT()
	{}

	virtual MREncodingDataType get_k_data( void );

	virtual void SampleFourierSpace( MREncodingDataType &i_data ) = 0;

protected:

	MREncodingDataType k_data_;


	
};


class FullySampledCartesianFFT: public aCartesianReadoutFFT{

public:
	FullySampledCartesianFFT();
	void SampleFourierSpace( MREncodingDataType &i_data);

};


class RadialPhaseEncodingFFT: public aCartesianReadoutFFT{

public:

	RadialPhaseEncodingFFT() {};

	void set_trajectory(TrajVessel &traj);
	void SampleFourierSpace( MREncodingDataType &i_data );

private:

	RPETrajectoryPreparation traj_prep_;
	
};

}