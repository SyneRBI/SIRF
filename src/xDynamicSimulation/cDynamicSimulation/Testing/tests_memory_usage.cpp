/* ================================================

Author: Johannes Mayer
Date: 2018.11.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include <time.h> 
#include <omp.h>

#include <ismrmrd/ismrmrd.h>

#include "auxiliary_input_output.h"

#include "tests_memory_usage.h"

using namespace std;
using namespace sirf;


void wait_for_time ( int const wait_time_s)
{

	clock_t endwait;
	endwait = clock () + wait_time_s * CLOCKS_PER_SEC ;
	while (clock() < endwait) {}
}


bool tests_memory::test_acquisition_memory( void )
{

	try
	{	
		AcquisitionsVector all_acquis;
		all_acquis.read( std::string( ISMRMRD_H5_TEST_PATH ));

		size_t num_iterations = 1000;

		for(size_t i=0; i<num_iterations; i++)
		{
			std::cout << "loopindex i " << i <<std::endl;
			AcquisitionsVector temp_dummy_vector;
			temp_dummy_vector.copy_acquisitions_info(all_acquis);
			

			for(size_t i_acq=0; i_acq<all_acquis.number(); i_acq++)
			{	
				ISMRMRD::Acquisition acq;
			    all_acquis.get_acquisition(i_acq, acq);
			    temp_dummy_vector.append_acquisition(acq);
			}

			ISMRMRD::Acquisition acq;	
			for(size_t j=0; j<temp_dummy_vector.items(); j++)
			{
				temp_dummy_vector.get_acquisition(j, acq);
			    std::cout << acq.scan_counter() << std::endl;
			}

		}

		return true;

	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}

bool tests_memory::test_acquisition_vector_memory( void )
{
	try
	{
		AcquisitionsVector all_acquis;// = mr_io::read_ismrmrd_acquisitions( ISMRMRD_H5_TEST_PATH );
		all_acquis.read( std::string( ISMRMRD_H5_TEST_PATH ));

		size_t num_acq = all_acquis.items();
		size_t num_iterations = 1000;

		std::vector< AcquisitionsVector > acq_vec_vector;

		for(size_t i=0; i<num_iterations; i++)
		{
			std::cout << " NUM IT "<<i<<std::endl;
			
			AcquisitionsVector temp_vec;
			temp_vec.copy_acquisitions_info( all_acquis ) ;

			for(size_t i_acq=0; i_acq<num_acq; i_acq++)
			{	
				auto sptr_acq = all_acquis.get_sptr_acquisition(i_acq);
				temp_vec.append_sptr_acquisition(sptr_acq);
			}
			acq_vec_vector.push_back( temp_vec );
		}


		return true;
	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}


bool tests_memory::test_ndarray_memory_managment( void )
{
	try
	{
		size_t const vol_dim_1D = 700;
		std::vector< size_t > volume_dims{vol_dim_1D,vol_dim_1D,vol_dim_1D};

		std::vector< size_t > empty_dims{0};

		int const num_iter = 100;

		for( int i=0; i<num_iter; i++)
		{
			std::cout << "Iterations # " << i <<std::endl;
			ISMRMRD::NDArray< double > temp_arr(volume_dims);
			
			#pragma omp parallel
			for( size_t i=0; i<temp_arr.getNumberOfElements(); i++)
			    *(temp_arr.begin() + i ) = 0.0;

			wait_for_time( 5 );

			temp_arr.resize( empty_dims );

		}
				    



		return true;
	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}