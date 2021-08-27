/* ================================================

Author: Johannes Mayer
Date: 2018.08.09
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include <math.h>
#include <time.h>
#include <random>
#include <stdexcept>
#include <stdlib.h>

#include "sirf/cDynamicSimulation/dynsim_noisegenerator.h"


using namespace sirf;

unsigned int aNoiseGenerator::generate_pseudo_seed( void )
{
	srand( time(NULL) );
	unsigned int seed = 0;
	while(seed == 0)
		seed = rand();

	return seed;
}

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


void GaussianNoiseGenerator::add_noise( MRAcquisitionData& ad ) 
{
	if( this->SNR_ > 0)
	{
		this->noise_width_kspace_ = noise_width_from_snr( ad );
		this->add_noise_to_data(ad);
	}
	else
	{
		this->add_noise_to_data( ad );
	}

}

float GaussianNoiseGenerator::noise_width_from_snr( MRAcquisitionData& ad )
{
	this->noise_width_img_ = this->signal_img_ / this->SNR_;

	size_t num_acquistions = ad.number();		
	
	if( num_acquistions <= 0)
		return 0.f;

	float const suggested_noise_width = this->noise_width_img_;

	return suggested_noise_width;
}

void GaussianNoiseGenerator::add_noise_to_data( MRAcquisitionData& ad ) 
{
	std::cout << "Adding gaussian noise of width " << this->noise_width_kspace_ <<std::endl;

	std::default_random_engine generator;
	generator.seed(this->generate_pseudo_seed());


	std::normal_distribution< float > gaussian_distribution( this->mean_noise_ , this->noise_width_kspace_ );

	for( size_t i_acq=0; i_acq<ad.number(); i_acq++)
	{
		ISMRMRD::Acquisition acq;
		ad.get_acquisition( i_acq, acq);

		for(size_t i_data_point=0; i_data_point<acq.getNumberOfDataElements(); i_data_point++)
		{
			*(acq.getDataPtr() + i_data_point) +=   complex_float_t(gaussian_distribution(generator), gaussian_distribution(generator)) / this->sequence_specific_scaling_;	
		}
	}

	std::cout << "finished." <<std::endl; ;
}
