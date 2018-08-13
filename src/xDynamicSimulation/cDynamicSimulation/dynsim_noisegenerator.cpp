/* ================================================

Author: Johannes Mayer
Date: 2018.08.09
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


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


void GaussianNoiseGenerator::add_noise( MRAcquisitionData& noisy_acquisition_data, MRAcquisitionData& noise_free_acquisition_data ) 
{

	std::default_random_engine generator;
	generator.seed( this->random_seed_ ); 
	std::normal_distribution< float > gaussian_distribution( this->mean_noise_ , this->width_noise_ );

	if(noisy_acquisition_data.number() != 0)
		throw std::runtime_error("Please give an empty container as noisy acquisitions. It will be filled with the noised data."); 

	ISMRMRD::Acquisition source_acq;

	for( size_t i_acq=0; i_acq<noise_free_acquisition_data.number(); i_acq++)
	{
		noise_free_acquisition_data.get_acquisition( i_acq, source_acq);
		ISMRMRD::Acquisition noisy_acq(source_acq);

		complex_float_t* const ptr_data =  noisy_acq.getDataPtr();
		
		for(size_t i_data_point=0; i_data_point<noisy_acq.getNumberOfDataElements(); i_data_point++)
		{
			*(ptr_data + i_data_point) +=  complex_float_t(gaussian_distribution(generator), gaussian_distribution(generator));	
		}

		noisy_acquisition_data.append_acquisition(noisy_acq);
	}
}