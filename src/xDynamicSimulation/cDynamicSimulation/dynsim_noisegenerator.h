/* ================================================

Author: Johannes Mayer
Date: 2018.08.09
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once


#include <memory>

#include "stir_types.h"

#include "stir_data_containers.h"
#include "gadgetron_data_containers.h"



typedef unsigned int SeedType;


class aNoiseGenerator{

public:

	aNoiseGenerator(){};
	aNoiseGenerator( SeedType const random_seed ): random_seed_(random_seed){};

	void const set_random_seed( SeedType const seed) {this->random_seed_ = seed;};
	SeedType const get_random_seed( void ) {return this->random_seed_;};

protected:
	SeedType random_seed_ = 0;

};



class PoissonNoiseGenerator: public aNoiseGenerator{

public:

	PoissonNoiseGenerator():aNoiseGenerator(){};

	std::shared_ptr<sirf::PETAcquisitionData> add_noise(sirf::PETAcquisitionData& acq_data){};

};

class GaussianNoiseGenerator: public aNoiseGenerator{

public:

	GaussianNoiseGenerator():aNoiseGenerator(){};
	GaussianNoiseGenerator(float const width_noise): aNoiseGenerator(), width_noise_(width_noise){};	

	std::shared_ptr<sirf::MRAcquisitionData> add_noise(sirf::MRAcquisitionData& acq_data){};

private:
	float width_noise_;

};