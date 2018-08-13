/* ================================================

Author: Johannes Mayer
Date: 2018.08.09
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "auxiliary_testing_functions.h"




#include "tests_noisegenerator.h"

using namespace sirf;


bool test_noisegen::test_add_poisson_noise( void )
{
	try
	{
			
		PETAcquisitionDataInFile noise_free_acq(FILENAME_DYNSIM_PET);
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
		AcquisitionsVector noise_free_acquisitions =  mr_io::read_ismrmrd_acquisitions( FILENAME_DYNSIM );
		AcquisitionsVector noisy_acquisitions;
		noisy_acquisitions.copy_acquisitions_info(noise_free_acquisitions);

		float const noise_width = 1;

		GaussianNoiseGenerator noise_gen( noise_width );
		noise_gen.add_noise(noisy_acquisitions, noise_free_acquisitions);

		std::string const filename_test_output = std::string(SHARED_FOLDER_PATH) + "test_gaussian_noise_generator.h5";		
		noisy_acquisitions.write(filename_test_output.c_str());

		return true;

	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}
}
