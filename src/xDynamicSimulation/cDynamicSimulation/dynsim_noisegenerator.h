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

	GaussianNoiseGenerator():aNoiseGenerator(){};
	GaussianNoiseGenerator(float const width_noise): aNoiseGenerator(), width_noise_(width_noise){};	

	void add_noise(sirf::MRAcquisitionData& noisy_acquisition_data, sirf::MRAcquisitionData& noise_free_acquisition_data );

private:
	static float constexpr mean_noise_ = 0.f;
	float width_noise_;

};