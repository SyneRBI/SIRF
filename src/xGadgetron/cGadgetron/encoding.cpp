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




void RPETrajectoryPreparation::set_and_check_trajectory( ISMRMRD::NDArray< TrajectoryPrecision > trajectory)
{

	std::vector<size_t> traj_dims = data_dims_from_ndarray< TrajectoryPrecision > ( trajectory );

	if( traj_dims[2] != 2)
		throw std::runtime_error("The trajectory must be formatted as a Ni x Nj x 2 array. Please give correct format.");
	
	this->traj_dims_.push_back( traj_dims[0] );
	this->traj_dims_.push_back( traj_dims[1] );

	this->traj_.create( this->traj_dims_ );
			
	size_t const num_traj_points = this->traj_dims_[0] * this->traj_dims_[1];

	for( size_t nr=0; nr<traj_dims_[0]; nr++)
		for( size_t na=0; na<traj_dims_[1]; na++)
		{
			TrajectoryPrecision traj_x = trajectory(nr, na, 0);
			TrajectoryPrecision traj_y = trajectory(nr, na, 1);

			size_t lin_index = nr*traj_dims_[1] + na;
			*(this->traj_.begin() + lin_index) = TrajectoryType2D(traj_x, traj_y);
		}	
}



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



void RadialPhaseEncodingFFT::set_trajectory(TrajectoryContainer &traj)
{
	this->traj_prep_.set_and_check_trajectory( traj );	
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

	std::vector<size_t> img_dims(data_dims.begin(),data_dims.begin() + 3) ;
    ho3DArray< complex_float_t > data_to_be_fftd( &img_dims );
		
	for( size_t i=0; i<num_elements; i++)
	{
		 auto val= *(i_data.begin() + i);
		 *(data_to_be_fftd.begin() + i) = val;
	}

	hoNDFFT< float >::instance()->fft1c( data_to_be_fftd );


	size_t num_slices = data_dims[0];
	auto traj_dims = this->traj_prep_.get_traj_dims();

	std::vector<size_t> output_data_size;
	output_data_size.push_back(num_slices);
	output_data_size.push_back(traj_dims[0]);
	output_data_size.push_back(traj_dims[1]);
	
	this->k_data_.resize(output_data_size);
	

	std::vector<size_t> slice_dims( data_dims.begin()+1, data_dims.begin()+3 ); 

	for(size_t i_slice=0; i_slice< num_slices; i_slice++)
	{
		ho2DArray< complex_float_t > sub_slice( &slice_dims );
		for(size_t y=0; y<data_dims[1]; y++)
		{
			for(size_t z=0; z<data_dims[2]; z++)
			{
				sub_slice(y, z) = data_to_be_fftd(i_slice, y, z);
			}
		}


		float const oversampling_factor = 1.f;
		float const kernel_size = 5.5f;


		hoNFFT_plan<float, 2> nufft_operator( from_std_vector<size_t, 2>(slice_dims) , oversampling_factor, kernel_size);
	
		hoNDArray< TrajectoryType2D > trajectory = this->traj_prep_.get_formatted_trajectory();
		
		nufft_operator.preprocess( trajectory );

		auto result = this->traj_prep_.get_formatted_output_container< complex_float_t >();
		
		auto identitiy_DCF = this->traj_prep_.get_formatted_identity_dcf< float >();

		nufft_operator.compute( sub_slice, result, identitiy_DCF, hoNFFT_plan<float, 2>::NFFT_FORWARDS_C2NC );


		for(size_t nr=0; nr<traj_dims[0]; nr++)
		{
			for(size_t na=0; na<traj_dims[1]; na++)
			{
				this->k_data_(i_slice, nr, na) = result(nr, na);
			}
		}

	}
}



