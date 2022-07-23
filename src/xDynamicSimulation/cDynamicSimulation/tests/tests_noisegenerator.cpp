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
	std::cout << "--- Running " << __FUNCTION__ << std::endl;
	try
	{
		std::stringstream fpath_pet_exampledata;
		fpath_pet_exampledata << SHARED_FOLDER_PATH << TESTDATA_PREFIX << "TemplateData/PET/simulated_data.hs"; 				

		PETAcquisitionDataInMemory noise_free_acq(fpath_pet_exampledata.str().c_str());
		PETAcquisitionDataInMemory noisy_acq(noise_free_acq);

		PoissonNoiseGenerator png;
		png.add_noise( noisy_acq, noise_free_acq);

		std::stringstream name_stream;
		name_stream << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__;

		std::cout << "Writing " << name_stream.str() << std::endl;
		noisy_acq.write(name_stream.str());

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
	std::cout << "--- Running " << __FUNCTION__ << std::endl;

	try
	{
		AcquisitionsVector av;
		av.read(ISMRMRD_H5_TEST_PATH);

		for(size_t i=0; i<av.number(); i++)
		{
			ISMRMRD::Acquisition acq;
			av.get_acquisition( i, acq );
			acq.resize( acq.number_of_samples() );

			auto it = acq.data_begin ();
			while( it != acq.data_end() )
			{
				*it = complex_float_t( 0.f, 0.f );	
				it++;
			} 
			av.set_acquisition(i, acq); 
		}

		float const SNR = 13;
		float const signal = 5;
		float const rpe_noise_scaling = 1.71558;


		GaussianNoiseGenerator noise_gen( SNR );
		noise_gen.set_signal_img(signal);
		noise_gen.set_sampling_specific_scaling( rpe_noise_scaling );

		noise_gen.add_noise(av);

		std::stringstream name_stream;
		name_stream << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__;
		av.write(name_stream.str());

		return true;

	}
	catch( std::runtime_error const &e)
	{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
	}
}
