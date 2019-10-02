/* ================================================

Author: Johannes Mayer
Date: 2018.11.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include <time.h> 
//TODO Johannes #include <omp.h>

#include <ismrmrd/ismrmrd.h>

#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "tests_memory_usage.h"

using namespace std;
using namespace sirf;


void wait_for_time ( int const wait_time_s)
{

	clock_t endwait;
	endwait = clock () + wait_time_s * CLOCKS_PER_SEC ;
	while (clock() < endwait) {}
}



bool tests_memory::test_sirf_free_acquisition( void )
{

	try
	{
		size_t num_iterations = 10000;

		size_t num_acquis = 256*128;
		uint16_t num_samples = 192;
		uint16_t num_channels = 32;
		
		ISMRMRD::Acquisition acq(num_samples, num_channels);

		bool construct_inside_loop = true;

		if( construct_inside_loop )
		{
			for(size_t i=0; i<num_iterations; i++)
			{
				std::cout << "loopindex i " << i <<std::endl;
		
				for(size_t i_acq=0; i_acq<num_acquis; i_acq++)
				{	
					ISMRMRD::Acquisition temp_acq;
				    temp_acq = acq;
				}
			}
		} 
		else if(!construct_inside_loop)
		{
			for(size_t i=0; i<num_iterations; i++)
			{
				std::cout << "loopindex i " << i <<std::endl;

				ISMRMRD::Acquisition temp_acq;
				for(size_t i_acq=0; i_acq<num_acquis; i_acq++)
				{
				    temp_acq = acq;
				} 
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


bool tests_memory::test_acquisition_memory( void )
{
	try
	{	
		AcquisitionsVector all_acquis;
		// all_acquis.read( std::string( ISMRMRD_H5_TEST_PATH ));
		all_acquis.read( std::string( ISMRMRD_H5_TEST_PATH ), false);

		size_t num_iterations = 1000;

		bool construct_inside_loop = false;

		
		if( construct_inside_loop )
		{
			for(size_t i=0; i<num_iterations; i++)
			{
				std::cout << "loopindex i " << i <<std::endl;
		
				for(size_t i_acq=0; i_acq<all_acquis.number(); i_acq++)
				{	
					ISMRMRD::Acquisition acq;
				    all_acquis.get_acquisition(i_acq, acq);
				}
			}
		} 
		else
		{
			for(size_t i=0; i<num_iterations; i++)
			{
				std::cout << "loopindex i " << i <<std::endl;

				ISMRMRD::Acquisition acq;
				for(size_t i_acq=0; i_acq<all_acquis.number(); i_acq++)
				{	
				    all_acquis.get_acquisition(i_acq, acq);
				}
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


bool tests_memory::test_downsizing_acquisition_memory( void )
{
	try
	{	
		AcquisitionsVector all_acquis;
		all_acquis.read( std::string( ISMRMRD_H5_TEST_PATH ));

		size_t num_iterations = 1000;

		AcquisitionsVector downsized_acquisitionsvector;
		downsized_acquisitionsvector.copy_acquisitions_info( all_acquis );
		for(size_t i_acq=0; i_acq<all_acquis.number(); i_acq++)
		{	
			ISMRMRD::Acquisition acq;
		    all_acquis.get_acquisition(i_acq, acq);
		    acq.resize(1,1);
			downsized_acquisitionsvector.append_acquisition(acq);
		}

		for(size_t i=0; i<num_iterations; i++)
		{
			std::cout << "loopindex i " << i <<std::endl;
			AcquisitionsVector temp_dummy_vector;
			temp_dummy_vector.copy_acquisitions_info(all_acquis);
		
			ISMRMRD::Acquisition acq;
			for(size_t i_acq=0; i_acq<downsized_acquisitionsvector.number(); i_acq++)
			{	
			    downsized_acquisitionsvector.get_acquisition(i_acq, acq);
				temp_dummy_vector.append_acquisition(acq);
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



bool tests_memory::test_acquisition_vector_ordering_memory( void )
{
	try
	{	
		AcquisitionsVector all_acquis;
		all_acquis.read( std::string( ISMRMRD_H5_TEST_PATH ));

		size_t num_iterations = 10000;

		for(size_t i=0; i<num_iterations; i++)
		{
			std::cout << "loopindex i " << i <<std::endl;
			std::cout << "Ordering " <<std::endl;
			
			all_acquis.order();
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
				auto sptr_acq = all_acquis.get_acquisition_sptr(i_acq);
				temp_vec.append_acquisition_sptr(sptr_acq);
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

bool tests_memory::tests_resizing_acquisition_memory( void )
{

	try
	{

		sirf::AcquisitionsVector av;
		av.read( std::string(ISMRMRD_H5_TEST_PATH)) ;

		size_t num_iterations = 10;

		for(size_t iteration=0; iteration<num_iterations; iteration++)
		{

			std::cout << "Waiting" << std::endl;
			wait_for_time( 5 );
			
			std::cout << " Resizing " <<std::endl;

			for(size_t i=0; i<av.number(); i++)
			{
				auto sptr_acq = av.get_acquisition_sptr(i);
				sptr_acq->resize( sptr_acq->number_of_samples(), (uint16_t)0, (uint16_t)0 );
				sptr_acq->resize( sptr_acq->number_of_samples(), (uint16_t)25, (uint16_t)0 );
				sptr_acq->resize( sptr_acq->number_of_samples(), (uint16_t)0, (uint16_t)0 );
			}

			auto some_sptr_acq = av.get_acquisition_sptr(0);
			std::cout << "Number of samples is now: " << some_sptr_acq->number_of_samples() <<std::endl;
			std::cout << "Number of active channels is now: " << some_sptr_acq->active_channels() <<std::endl;
			std::cout << "Data size is now: " << some_sptr_acq->getDataSize() << std::endl;

			std::cout << "Waiting again" <<std::endl;
			wait_for_time( 5 );
			

			std::cout << "Resizing again" <<std::endl;
			for(size_t i=0; i<av.number(); i++)
			{
				auto sptr_acq = av.get_acquisition_sptr(i);
				sptr_acq->resize( sptr_acq->number_of_samples(), 4 );

			}

			wait_for_time( 5 );
			std::cout << "Waiting again" <<std::endl;
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


bool  tests_memory::tests_VD_h5_file_content( void )
{

try
	{

		sirf::AcquisitionsVector av;
		av.read( std::string(ISMRMRD_H5_TEST_PATH)) ;

		
		std::cout << "# of acquis in vector " << av.number() << std::endl;

		return true;
	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}