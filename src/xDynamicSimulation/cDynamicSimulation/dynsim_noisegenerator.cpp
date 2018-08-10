/* ================================================

Author: Johannes Mayer
Date: 2018.08.09
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



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