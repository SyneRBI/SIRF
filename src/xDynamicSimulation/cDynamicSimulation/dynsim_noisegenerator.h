/* ================================================

Author: Johannes Mayer
Date: 2018.08.09
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once


#include <memory>

#include <stir/GeneralisedPoissonNoiseGenerator.h>

#include "stir_types.h"

#include "stir_data_containers.h"
#include "gadgetron_data_containers.h"

#define RPE_NOISE_SCALING 1.71558f


typedef unsigned int SeedType;


class aNoiseGenerator{

public:

	aNoiseGenerator(){};
	aNoiseGenerator( SeedType const random_seed ): random_seed_(random_seed){}
	virtual void set_random_seed( SeedType const seed) {this->random_seed_ = seed;};
	SeedType const get_random_seed( void ) {return this->random_seed_;};

protected:
	SeedType random_seed_ = 1;

};



class PoissonNoiseGenerator: public aNoiseGenerator{

public:

	PoissonNoiseGenerator():aNoiseGenerator(), stir_noise_gen_(1.0f,true)
	{
		this->stir_noise_gen_.seed(this->random_seed_);
	}


	virtual void set_random_seed( SeedType const seed) 
	{
		this->random_seed_ = seed;
		this->stir_noise_gen_.seed(this->random_seed_);
	}

	void add_noise( sirf::PETAcquisitionData& noisy_acq, sirf::PETAcquisitionData& noise_free_acq);

private:
	stir::GeneralisedPoissonNoiseGenerator stir_noise_gen_;	

};

class GaussianNoiseGenerator: public aNoiseGenerator{

public:

	GaussianNoiseGenerator():aNoiseGenerator(), noise_width_img_(0.f){};
	GaussianNoiseGenerator(float const SNR): aNoiseGenerator(), SNR_(SNR){};	

	void set_signal_img(float const signal){ this->signal_img_ = signal; };
	void set_SNR(float const SNR){ this->SNR_ = SNR; };
	void set_sampling_specific_scaling( float const scaling) { this->sequence_specific_scaling_ = scaling;};

	void add_noise( sirf::AcquisitionsVector& acquisition_vector );
	

private:
	static float constexpr mean_noise_ = 0.f;
	
	float SNR_ = -1;
	float signal_img_ = 0.f;
	float noise_width_img_ =0.f;
	float noise_width_kspace_ =0.f;

	float sequence_specific_scaling_ = 1.f;	
	
	void add_noise_to_data( sirf::AcquisitionsVector& acquisition_vector );
	float noise_width_from_snr( sirf::AcquisitionsVector& acquisition_vector );

};