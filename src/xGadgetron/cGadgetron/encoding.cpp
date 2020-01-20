/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include <stdexcept>
#include <algorithm>

#include "sirf/cGadgetron/encoding.h"


#include <gadgetron/hoNDArray.h>
#include <gadgetron/ho2DArray.h>
#include <gadgetron/ho3DArray.h>

#include <gadgetron/hoNDFFT.h>
#include <gadgetron/hoNFFT.h>


#include <gadgetron/hoNDArray_utils.h>

#include <gadgetron/vector_td_utilities.h>
#include <gadgetron/vector_td_io.h>


using namespace sirf;
using namespace Gadgetron;
using namespace ISMRMRD;
using ISMRMRD::ISMRMRD_NDARRAY_MAXDIM;





void aTrajectoryContainer::overwrite_ismrmrd_trajectory_info(ISMRMRD::IsmrmrdHeader& hdr) 
{
	std::vector<ISMRMRD::Encoding> current_encoding = hdr.encoding;
	std::vector<ISMRMRD::Encoding> overwritten_encoding;
	
	for( int i=0; i<current_encoding.size(); i++)
	{
		ISMRMRD::Encoding enc = current_encoding[i];
		enc.trajectory = traj_type_;
		overwritten_encoding.push_back(enc);
	}

	hdr.encoding = overwritten_encoding;
}


void aTrajectoryContainer::overwrite_ismrmrd_trajectory_info(std::string& serialized_header) 
{
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(serialized_header.c_str() , hdr);

	this->overwrite_ismrmrd_trajectory_info( hdr );

	std::stringstream updated_serialized_header_stream;
    ISMRMRD::serialize( hdr, updated_serialized_header_stream);
    serialized_header = updated_serialized_header_stream.str();
}

void aTrajectoryContainer::set_header(ISMRMRD::IsmrmrdHeader hdr)
{
	this->hdr_ = hdr;
	this->overwrite_ismrmrd_trajectory_info(this->hdr_);
}

void aTrajectoryContainer::set_trajectory( TrajVessel trajectory )
{
	this->traj_ = trajectory;
}

TrajVessel aTrajectoryContainer::get_trajectory( void )
{
	return this->traj_;
}

ISMRMRD::TrajectoryType aTrajectoryContainer::get_traj_type( void )
{
    // throw std::runtime_error("aTrajectoryContainer::get_traj_type: This should return ISMRMRD::TrajectoryType, not std::string");
    // TODO Johannes Uncomment below when fixed.
    return this->traj_type_;
}

void aTrajectoryContainer::norm_trajectory( void )
{
	TrajPrecision traj_max = this->get_traj_max_abs();
	
	for(int j=0; j<this->traj_.getNumberOfElements(); j++)
		*(this->traj_.begin() + j) /= (2*traj_max); 

}



TrajPrecision RPETrajectoryContainer::get_traj_max_abs( void )
{
	auto traj_dims = this->traj_.getDims();	

	TrajPrecision maximum_traj_value = 0;

	for(int na=0; na<traj_dims[1]; na++)
	for(int nr=0; nr<traj_dims[0]; nr++)
	{

		std::complex< TrajPrecision > x( this->traj_(nr,na,0), this->traj_(nr,na,1) );
		TrajPrecision abs_x = std::abs(x);
		if(abs_x  > maximum_traj_value)
			maximum_traj_value = abs_x;
	}

	return maximum_traj_value;
}


void RPETrajectoryContainer::compute_trajectory()
{

	std::cout << "Computing trajectory" << std::endl;

	std::vector< Encoding > all_encodings = this->hdr_.encoding; 

	if( all_encodings.size() != 1 )
		throw LocalisedException("Your header file contains zero or more than one encodings. Please pass one with only one encoding.", __FILE__, __LINE__);

	Encoding enc = all_encodings[0];
	EncodingSpace enc_space = enc.encodedSpace;
	
	MatrixSize encoding_mat_size = enc_space.matrixSize;

	unsigned short NRadial = encoding_mat_size.y;
	unsigned short NAngles = encoding_mat_size.z;

	std::vector<size_t> traj_dims{NRadial, NAngles, 2}; 

   	this->traj_.resize(traj_dims);

	for( unsigned na=0; na<NAngles; na++){
	for( unsigned nr=0; nr<NRadial; nr++)
	{

		int const r_pos = nr - NRadial /2;
		float const ang_pos = na*M_PI/ NAngles;
			
		float const nx = r_pos * cos( ang_pos );
		float const ny = r_pos * sin( ang_pos );

		this->traj_(nr, na, 0) = nx;
		this->traj_(nr, na, 1) = ny;
	}}

	this->norm_trajectory();
}

void RPETrajectoryContainer::set_acquisition_trajectory(ISMRMRD::Acquisition& acq)
{
	if ( this->traj_.getNumberOfElements() == 0)
		throw LocalisedException("No trajectory has been set or computed. I suggest you execute compute_trajectory first in order to be", __FILE__, __LINE__);


	auto num_samples = acq.number_of_samples();
	auto active_channels = acq.active_channels();
	auto trajectory_dimensions = 3;

	acq.resize(num_samples, active_channels, trajectory_dimensions);

	uint16_t const enc_step_1 = acq.getHead().idx.kspace_encode_step_1;
	uint16_t const enc_step_2 = acq.getHead().idx.kspace_encode_step_2;
	
	TrajPrecision const ky = this->traj_( enc_step_1, enc_step_2, 0);
	TrajPrecision const kz = this->traj_( enc_step_1, enc_step_2 ,1); 

	for( unsigned i=0; i<num_samples; i++)
	{
		float const readout_traj = (-(float)num_samples/2.f + (float)i ) / (float)num_samples;
		
		acq.traj(0, i )		= readout_traj;
		acq.traj(1, i ) 	= ky;
		acq.traj(2, i ) 	= kz;

	}
}


void RPEInterleavedTrajectoryContainer::compute_trajectory()
{
	std::cout << "Computing trajectory" << std::endl;

	std::vector< Encoding > all_encodings = this->hdr_.encoding; 

	if( all_encodings.size() != 1 )
		throw LocalisedException("Your header file contains zero or more than one encodings. Please pass one with only one encoding.", __FILE__, __LINE__);

	Encoding enc = all_encodings[0];
	EncodingSpace enc_space = enc.encodedSpace;
	
	MatrixSize encoding_mat_size = enc_space.matrixSize;

	unsigned short NRadial = encoding_mat_size.y;
	unsigned short NAngles = encoding_mat_size.z;

	std::vector<size_t> traj_dims{NRadial, NAngles, 2}; 

   	this->traj_.resize(traj_dims);

   	std::vector<float> radial_shift{0.f, 2.f, 1.f, 3.f};

	
	for( unsigned na=0; na<NAngles; na++){
	for( unsigned nr=0; nr<NRadial; nr++)
	{

		float r_pos = (float)nr - (float)NRadial/2.f; 
		r_pos = (r_pos == 0) ? r_pos : r_pos + 0.25f * radial_shift[ na % 4 ];
		float const ang_pos = na*M_PI/ NAngles;
				
		float const nx = r_pos * cos( ang_pos );
		float const ny = r_pos * sin( ang_pos );

		this->traj_(nr, na, 0) = nx;
		this->traj_(nr, na, 1) = ny;
	}}

	this->norm_trajectory();
}


#define GOLDENRATIO (1 + sqrt(5.f))/2.f
// #define GOLDENANGLE M_PI*(3 - sqrt(5))
#define GOLDENANGLE M_PI*(GOLDENRATIO-1)


void RPEInterleavedGoldenCutTrajectoryContainer::compute_trajectory()
{
	
	std::cout << "Computing trajectory" << std::endl;

	std::vector< Encoding > all_encodings = this->hdr_.encoding; 

	if( all_encodings.size() != 1 )
		throw LocalisedException("Your header file contains zero or more than one encodings. Please pass one with only one encoding.", __FILE__, __LINE__);

	Encoding enc = all_encodings[0];
	EncodingSpace enc_space = enc.encodedSpace;

	MatrixSize encoding_mat_size = enc_space.matrixSize;

	unsigned short NRadial = encoding_mat_size.y;
	unsigned short NAngles = encoding_mat_size.z;

	std::vector<size_t> traj_dims{NRadial, NAngles, 2}; 

   	this->traj_.resize(traj_dims);

   	std::vector<float> radial_shift{0.f, 2.f, 1.f, 3.f};

	for( unsigned na=0; na<NAngles; na++)
	for( unsigned nr=0; nr<NRadial; nr++)
	{

		float r_pos = (float)nr - (float)NRadial/2.f; 
		r_pos = (r_pos == 0) ? r_pos : r_pos + 0.25f * radial_shift[ na % 4 ];
		float const ang_pos = na*GOLDENANGLE;
				
		float const nx = r_pos * cos( ang_pos );
		float const ny = r_pos * sin( ang_pos );

		this->traj_(nr, na, 0) = nx;
		this->traj_(nr, na, 1) = ny;
	}
	this->norm_trajectory();
}


void RPESuperInterleavedGoldenCutTrajectoryContainer::compute_trajectory()
{
	std::cout << "Computing trajectory" << std::endl;

	std::vector< Encoding > all_encodings = this->hdr_.encoding; 

	if( all_encodings.size() != 1 )
		throw LocalisedException("Your header file contains zero or more than one encodings. Please pass one with only one encoding.", __FILE__, __LINE__);

	Encoding enc = all_encodings[0];
	EncodingSpace enc_space = enc.encodedSpace;

	MatrixSize encoding_mat_size = enc_space.matrixSize;

	unsigned short NRadial = encoding_mat_size.y;
	unsigned short NAngles = encoding_mat_size.z;

	std::vector<size_t> traj_dims{NRadial, NAngles, 2}; 

   	this->traj_.resize(traj_dims);

	std::vector<PermutIdxContainer> sigma;

	for( unsigned na=0; na<NAngles; na++)
	{
		PermutIdxContainer angle_with_index;
		angle_with_index.first = fmod(na * GOLDENANGLE, M_PI);
		angle_with_index.second = na;

		sigma.push_back(angle_with_index);
	}

	std::sort(sigma.begin(), sigma.end(), this->comparefun);

	std::vector<size_t> gc_sorted_lut(NAngles);
	for( unsigned na=0; na<NAngles; na++)
		gc_sorted_lut[ sigma[na].second ] = na;

	float const angle_increment = GOLDENANGLE;

	for( unsigned na=0; na<NAngles; na++)
	{
		float const ang_pos =  fmod(na*angle_increment,M_PI); 
		float const curr_shift = 0.5 * (2*fmod((gc_sorted_lut[na])*GOLDENRATIO, 1) -1);
		
		for( unsigned nr=0; nr<NRadial; nr++)
		{

			float r_pos = (float)nr - (float)NRadial/2.f; 
			r_pos = (nr == 0) ? 0 : r_pos + curr_shift;
					
			float const nx = r_pos * cos( ang_pos );
			float const ny = r_pos * sin( ang_pos );

			this->traj_(nr, na, 0) = nx;
			this->traj_(nr, na, 1) = ny;
		}
	}

	this->norm_trajectory();
}



void RPETrajectoryPreparation::set_and_check_trajectory( TrajVessel& trajectory)
{
	
	this->lut_idx_to_traj_.clear();

	std::vector<size_t> traj_dims = data_dims_from_ndarray< TrajPrecision > ( trajectory );

	if( traj_dims[2] != 2)
		throw std::runtime_error("The trajectory must be formatted as a Ni x Nj x 2 array. Please give correct format.");
	
	this->traj_dims_.push_back( traj_dims[0] );
	this->traj_dims_.push_back( traj_dims[1] );

	this->full_traj_dims_ = traj_dims_;

	this->traj_.create( this->traj_dims_ );
			
	size_t const num_traj_points = this->traj_dims_[0] * this->traj_dims_[1];

	this->lut_idx_to_traj_.resize( num_traj_points );

	for( size_t na=0; na<traj_dims_[1]; na++)
	for( size_t nr=0; nr<traj_dims_[0]; nr++)
	{
		TrajPrecision traj_x = trajectory(nr, na, 0);
		TrajPrecision traj_y = trajectory(nr, na, 1);

		size_t const lin_index = na*traj_dims_[0] + nr;
		this->lut_idx_to_traj_[ lin_index ] = std::make_pair(nr, na);

		*(this->traj_.begin() + lin_index) = TrajectoryType2D(traj_x, traj_y);
	}	
}

void RPETrajectoryPreparation::set_kspace_subset( sirf::MRAcquisitionData &ad )
{
	if( this->traj_.get_number_of_elements() < 1 )
		throw std::runtime_error( "Please set trajectory before trying to get a kspace subset");

	size_t const num_traj_points = ad.number();

	this->lut_idx_to_traj_.clear();
	this->lut_idx_to_traj_.resize( num_traj_points );

	std::vector< size_t > subset_traj_dims{num_traj_points, 1};
	
	Gadgetron::hoNDArray< TrajectoryType2D > subset_traj;	
	subset_traj.create( subset_traj_dims );

	for( size_t i=0; i<num_traj_points; i++)
	{

		auto sptr_acq = ad.get_acquisition_sptr( i );

		size_t const nr = (size_t)sptr_acq->getHead().idx.kspace_encode_step_1;
		size_t const na = (size_t)sptr_acq->getHead().idx.kspace_encode_step_2;
		size_t const lin_index = na*this->traj_dims_[0] + nr;
		
		*(subset_traj.begin() + i) = *(this->traj_.begin() + lin_index);
		
		lut_idx_to_traj_[i] = std::make_pair(nr, na);
	}

	this->traj_ = subset_traj;
	this->traj_dims_ = subset_traj_dims;
}

std::pair< size_t, size_t >  RPETrajectoryPreparation::get_traj_index_pair( size_t const idx ) const
{
	return this->lut_idx_to_traj_[idx];	
}

MREncodingDataType aCartesianReadoutFFT::get_k_data( void )
{
	return this->k_data_;	
}


FullySampledCartesianFFT::FullySampledCartesianFFT():
aCartesianReadoutFFT()
{
}

void FullySampledCartesianFFT::SampleFourierSpace( MREncodingDataType &i_data)
{

	size_t const num_elements = i_data.getNumberOfElements();

	std::vector<size_t> data_dims = data_dims_from_ndarray<complex_float_t>(i_data);

	Gadgetron::hoNDArray< complex_float_t > data_to_be_fftd(data_dims);
	
	for( size_t i=0; i<num_elements; i++)
		*(data_to_be_fftd.begin() + i) = *(i_data.begin() + i);

	// transform first three dims
	hoNDFFT< float >::instance()->fft3c( data_to_be_fftd);
	
	this->k_data_ = MREncodingDataType(data_dims);

	for( size_t i=0; i<num_elements; i++)
	{
		*(this->k_data_.begin() + i) = *(data_to_be_fftd.begin() + i);
	}

}



void RadialPhaseEncodingFFT::set_trajectory(TrajVessel &traj)
{
	this->traj_prep_.set_and_check_trajectory( traj );	
}

void RadialPhaseEncodingFFT::set_kspace_subset( sirf::MRAcquisitionData& ad )
{
	this->traj_prep_.set_kspace_subset( ad );
}


void RadialPhaseEncodingFFT::SampleFourierSpace( MREncodingDataType &i_data)
{

	std::vector<size_t> data_dims = data_dims_from_ndarray< complex_float_t >(i_data);

	size_t num_slices = data_dims[0];
	size_t const num_coils = data_dims[3];
		
	hoNDArray< complex_float_t > data_to_be_fftd( data_dims );


	auto traj_dims = this->traj_prep_.get_full_traj_dims();

	size_t const num_elements = i_data.getNumberOfElements();
	
	#pragma omp parallel
	for( size_t i=0; i<num_elements; i++)
		*(data_to_be_fftd.begin() + i) = *(i_data.begin() + i);
	
	// memcpy( data_to_be_fftd.begin(), i_data.begin(), i_data.getDataSize());

	
	hoNDFFT< float >::instance()->fft1c( data_to_be_fftd );
	
	size_t const oversampling_factor = 2;
	size_t const kernel_size = 2;	//must keep integers! nfft instable with floats

	std::vector<size_t> slice_dims( data_dims.begin()+1, data_dims.begin()+3 ); 

	size_t const usos_factor = 4;

	size_t const Nr = traj_dims[0]*usos_factor ;
	
	std::vector<size_t> cropped_slice_dims {Nr, Nr}; 
	std::vector<size_t> crop_offset_idx{data_dims[1]/2 - Nr/2, data_dims[2]/2 - Nr/2};

	Gadgetron::uint64d2 cropped_slice_dimension = from_std_vector< size_t, 2>(cropped_slice_dims);
	Gadgetron::uint64d2 crop_offset = from_std_vector< size_t, 2>(crop_offset_idx);

	if ( Nr > data_dims[1] || Nr > data_dims[2] )
	{
		std::cout << "Nr = " << Nr << std::endl;
		std::cout << "data_dims[1,2] = " << data_dims[1] << "," << data_dims[2] << std::endl;
		throw std::runtime_error( "You passed a volume with less voxels than radial encoding points in your MR template file."); 
	}
	
	std::vector<size_t> oversampled_slice_dims;
	
	for (size_t i=0; i<cropped_slice_dims.size(); i++)
		oversampled_slice_dims.push_back( oversampling_factor * cropped_slice_dims[i] );	
		
	Gadgetron::uint64d2 gridder_img_dimensions = from_std_vector<size_t, 2>(oversampled_slice_dims);
	hoNFFT_plan<float, 2> nufft_operator( gridder_img_dimensions , (float)1, (float)kernel_size);// dont change the osf -> nfft instable wrt segfaults


	hoNDArray< TrajectoryType2D > trajectory = this->traj_prep_.get_formatted_trajectory();
	
	size_t const num_traj_elem = trajectory.get_number_of_elements();
	
	nufft_operator.preprocess( trajectory );
	bool found_bad_val = false;

	for(size_t i_coil=0; i_coil<num_coils; i_coil++)
	{
		for(size_t i_slice=0; i_slice< num_slices; i_slice++)
		{
			ho2DArray< complex_float_t > sub_slice( &slice_dims );

			for(size_t z=0; z<data_dims[2]; z++)
			{
				for(size_t y=0; y<data_dims[1]; y++)
				{
					size_t const lin_index_access_4D = ((i_coil * data_dims[2] + z)*data_dims[1] + y)*num_slices + i_slice;

					sub_slice(y, z) = data_to_be_fftd [ lin_index_access_4D ];		
				}
			}

			ho2DArray< complex_float_t > cropped_subslice( &cropped_slice_dims );
			Gadgetron::crop<complex_float_t, 2>(crop_offset, cropped_slice_dimension, &sub_slice, &cropped_subslice);

		    ho2DArray<complex_float_t> padded_sub_slice;
		    Gadgetron::pad<complex_float_t,  2>(gridder_img_dimensions, &cropped_subslice, &padded_sub_slice, true, 0.f);

			Gadgetron::hoNDArray< complex_float_t > result = this->traj_prep_.get_formatted_output_container< complex_float_t >();
			Gadgetron::hoNDArray< float > identitiy_DCF = this->traj_prep_.get_formatted_identity_dcf< float >();

			nufft_operator.compute( padded_sub_slice, result, identitiy_DCF, hoNFFT_plan<float, 2>::NFFT_FORWARDS_C2NC );

			std::vector<size_t> k_space_size;
			k_space_size.push_back(num_slices);
			k_space_size.push_back(traj_dims[0]);
			k_space_size.push_back(traj_dims[1]);
			k_space_size.push_back(num_coils);
			this->k_data_.resize(k_space_size);

			for(size_t i_elem=0; i_elem< result.get_number_of_elements(); i_elem++) 	
			{
				auto radius_angle_idx_pair = this->traj_prep_.get_traj_index_pair( i_elem );

				size_t const nr = radius_angle_idx_pair.first;
				size_t const na = radius_angle_idx_pair.second;

				this->k_data_(i_slice, nr, na, i_coil) = result[i_elem];
				
				if( std::abs(result[i_elem]) > 1e9 && found_bad_val == false)
				{
					std::cout << "Potentially large value: " << result[i_elem] <<std::endl;
					found_bad_val = true;
				}
			}
			
		}	
	}
}




