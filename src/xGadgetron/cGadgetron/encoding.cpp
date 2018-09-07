/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include <stdexcept>

#include "encoding.h"

#include <gadgetron/hoNDArray.h>
#include <gadgetron/ho2DArray.h>
#include <gadgetron/ho3DArray.h>

#include <gadgetron/hoNDFFT.h>
#include <gadgetron/hoNFFT.h>

#include <gadgetron/vector_td_utilities.h>

using namespace Gadgetron;
using ISMRMRD::ISMRMRD_NDARRAY_MAXDIM;


ISMRMRD::NDArray<complex_float_t> aCartesianReadoutFFT::get_k_data( void )
{
	return this->k_data_;	
}


FullySampledCartesianFFT::FullySampledCartesianFFT():
aCartesianReadoutFFT()
{
}

void FullySampledCartesianFFT::SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data)
{

	size_t const num_elements = i_data.getNumberOfElements();

	std::vector<size_t> data_dims = data_dims_from_ndarray<complex_float_t>(i_data);

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





void RadialPhaseEncodingFFT::set_trajectory(TrajectoryType &traj)
{
	this->traj_ = traj;
}


void RadialPhaseEncodingFFT::SampleFourierSpace( ISMRMRD::NDArray<complex_float_t> i_data)
{
	
	size_t const num_elements = i_data.getNumberOfElements();

	std::vector<size_t> data_dims = data_dims_from_ndarray< complex_float_t >(i_data);

	bool is_3d_data = true;
	for(int i=3; i<ISMRMRD_NDARRAY_MAXDIM; i++)
		is_3d_data *= (data_dims[i] == 1);
	
	if( ! is_3d_data )
		throw std::runtime_error("Please pass 3D data to the RPE fourier transform");

    ho3DArray< complex_float_t > data_to_be_fftd( &data_dims );
		
	for( size_t i=0; i<num_elements; i++)
	{
		 auto val= *(i_data.begin() + i);
		 *(data_to_be_fftd.begin() + i) = val;
	}

	hoNDFFT< float >::instance()->fft1c( data_to_be_fftd );


	std::vector<size_t> slice_dims( data_dims.begin()+1, data_dims.begin()+3 ); 

	for(size_t i_slice=0; i_slice<data_dims[0]; i_slice++)
	{
		ho2DArray< float_complext > sub_slice( &slice_dims );
		for(size_t y=0; y<data_dims[1]; y++)
		{
			for(size_t z=0; z<data_dims[2]; z++)
			{
				sub_slice(y, z) = data_to_be_fftd(i_slice, y, z);
			}
		}

		hoNFFT_plan<float, 2> plan( from_std_vector<size_t, 2>(slice_dims) , 2.f, 1.f);
		
		std::vector<size_t> spatial_traj_dims = data_dims_from_ndarray< std::vector<float> >( this->traj_ );

		ho2DArray< floatd2 > trajectory( &spatial_traj_dims );
		
		size_t const num_traj_points = this->traj_.getNumberOfElements();

		for( size_t i=0; i<num_traj_points; i++)
			*(trajectory.begin() + i) = from_std_vector<float, 2>( *(this->traj_.begin() + i) );		

		plan.preprocess( trajectory );

		ho2DArray<float_complext> result; result.create(spatial_traj_dims[0], spatial_traj_dims[1]);
		plan.compute( sub_slice, result,hoNFFT_plan<float, 2>::NFFT_FORWARD_C2NC );

	}
}
