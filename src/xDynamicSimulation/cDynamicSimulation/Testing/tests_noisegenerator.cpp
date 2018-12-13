/* ================================================

Author: Johannes Mayer
Date: 2018.08.09
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "auxiliary_testing_functions.h"

#include <sstream>


#include "tests_noisegenerator.h"

using namespace sirf;


bool test_noisegen::test_add_poisson_noise( void )
{
	try
	{
			
		PETAcquisitionDataInFile noise_free_acq(FILENAME_STATICSIM_PET);
		PETAcquisitionDataInMemory noisy_acq(noise_free_acq);

		PoissonNoiseGenerator png;
		png.add_noise( noisy_acq, noise_free_acq);

		std::string const filename_test_output = std::string(SHARED_FOLDER_PATH) + "test_poisson_noise_generator.hs";		

		std::cout << "Writing " << filename_test_output << std::endl;
		noisy_acq.write(filename_test_output.c_str());

		return true;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}


bool test_noisegen::test_add_gaussian_noise( void )
{
	try
	{
		AcquisitionsVector av;
		av.read(ISMRMRD_H5_TEST_PATH);

		for(size_t i=0; i<av.number(); i++)
		{
			auto sptr_ac = av.get_acquisition_sptr( i );
			sptr_ac->resize( sptr_ac->number_of_samples() );

			auto it = sptr_ac->data_begin ();
			while( it != sptr_ac->data_end() )
			{
				*it = complex_float_t( 0.f, 0.f );	
				it++;
			}  
		}

		float const SNR = 13;
		float const signal = 5;
		float const rpe_noise_scaling = 1.71558;


		GaussianNoiseGenerator noise_gen( SNR );
		noise_gen.set_signal_img(signal);
		noise_gen.set_sequence_specific_scaling( rpe_noise_scaling );

		noise_gen.add_noise(av);

		std::stringstream filename_test_output;
		filename_test_output << std::string(SHARED_FOLDER_PATH) << "test_gaussian_noise_generator_SNR_" << SNR <<"_signal_" << signal  << ".h5";		
		av.write(filename_test_output.str().c_str());

		return true;

	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}
}
