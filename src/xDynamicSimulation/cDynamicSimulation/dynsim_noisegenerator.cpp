/* ================================================

Author: Johannes Mayer
Date: 2018.08.09
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include <math.h>
#include <random>
#include <stdexcept>


#include "dynsim_noisegenerator.h"


using namespace sirf;


void PoissonNoiseGenerator::add_noise( PETAcquisitionData& noisy_acq, PETAcquisitionData& noise_free_acq)
{
	try
	{
		auto sptr_noise_free_data = noise_free_acq.data();
		auto sptr_noisy_data = noisy_acq.data();

		this->stir_noise_gen_.generate_random( *(sptr_noisy_data), *(sptr_noise_free_data) );	
		noisy_acq.set_data( sptr_noisy_data );
	}
	catch(...)
	{
		std::cout << "Something went wrong in adding poisson noise" << std::endl;
	}
}

void GaussianNoiseGenerator::add_noise( AcquisitionsVector& acquisition_vector ) 
{
	if( this->SNR_ > 0)
	{
		this->width_noise_ = noise_width_from_snr( acquisition_vector );
		this->add_noise_to_data(acquisition_vector);
	}
	else
	{
		this->add_noise_to_data( acquisition_vector );
	}

}

float GaussianNoiseGenerator::noise_width_from_snr( AcquisitionsVector& acquisition_vector )
{
	size_t number_sampled_k_space_points = 0;
	float sum_max_signals = 0;

	size_t num_acquistions = acquisition_vector.number();

	for( size_t i_acq=0; i_acq<num_acquistions; i_acq++)
	{
		ISMRMRD::Acquisition acq;
		acquisition_vector.get_acquisition( i_acq, acq);
		
		size_t const num_data_points_curr_acquisition = acq.getNumberOfDataElements();
		
		complex_float_t* const ptr_data =  acq.getDataPtr();
		

		float curr_max_signal = 0;

		for(size_t i_data_point=0; i_data_point<num_data_points_curr_acquisition; i_data_point++)
		{
			curr_max_signal = ( std::abs( *(ptr_data + i_data_point) ) > curr_max_signal)? std::abs( *(ptr_data + i_data_point) ) : curr_max_signal;
		}
		sum_max_signals += curr_max_signal;


	}


	float const suggested_noise_width = sum_max_signals / ( (float)num_acquistions * sqrt(2.f) * this->SNR_ );

	std::cout << "Noise width: " <<suggested_noise_width << std::endl;
	return suggested_noise_width;
}

void GaussianNoiseGenerator::add_noise_to_data( AcquisitionsVector& acquisition_vector ) 
{

	std::cout << "Adding noise ..." ;
	acquisition_vector.set_as_template();

	std::default_random_engine generator;
	generator.seed( this->random_seed_ ); 
	std::normal_distribution< float > gaussian_distribution( this->mean_noise_ , this->width_noise_ );


	for( size_t i_acq=0; i_acq<acquisition_vector.number(); i_acq++)
	{
		auto sptr_acquis = acquisition_vector.get_sptr_acquisition( i_acq );

		for(size_t i_data_point=0; i_data_point<sptr_acquis->getNumberOfDataElements(); i_data_point++)
		{
			*(sptr_acquis->getDataPtr() + i_data_point) +=  complex_float_t(gaussian_distribution(generator), gaussian_distribution(generator));	
		}
	}
	std::cout << "finished." <<std::endl; ;
}

// void GaussianNoiseGenerator::add_noise( MRAcquisitionData& noisy_acquisition_data, MRAcquisitionData& noise_free_acquisition_data ) 
// {

// 	std::default_random_engine generator;
// 	generator.seed( this->random_seed_ ); 
// 	std::normal_distribution< float > gaussian_distribution( this->mean_noise_ , this->width_noise_ );

// 	if(noisy_acquisition_data.number() != 0)
// 		throw std::runtime_error("Please give an empty container as noisy acquisitions. It will be filled with the noised data."); 


// 	for( size_t i_acq=0; i_acq<noise_free_acquisition_data.number(); i_acq++)
// 	{
// 		ISMRMRD::Acquisition source_acq;
// 		noise_free_acquisition_data.get_acquisition( i_acq, source_acq);
		
// 		ISMRMRD::Acquisition noisy_acq(source_acq);

// 		complex_float_t* const ptr_data =  noisy_acq.getDataPtr();
		
// 		for(size_t i_data_point=0; i_data_point<noisy_acq.getNumberOfDataElements(); i_data_point++)
// 		{
// 			*(ptr_data + i_data_point) +=  complex_float_t(gaussian_distribution(generator), gaussian_distribution(generator));	
// 		}

// 		noisy_acquisition_data.append_acquisition(noisy_acq);
// 	}
// }