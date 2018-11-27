/* ================================================

Author: Johannes Mayer
Date: 2018.11.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "tests_memory_usage.h"

#include "auxiliary_input_output.h"

#include <ismrmrd/ismrmrd.h>

using namespace std;
using namespace sirf;



bool tests_memory::test_acquisition_memory( void )
{

	try
	{	
		AcquisitionsVector all_acquis = mr_io::read_ismrmrd_acquisitions( ISMRMRD_H5_TEST_PATH );

		PtrAcquisVect sptr_vec;

		for(size_t i_acq=0; i_acq<all_acquis.number(); i_acq++)
		{
			ISMRMRD::Acquisition acq;
		    all_acquis.get_acquisition(i_acq, acq);
		    sptr_vec.append_sptr_acq( std::make_shared< ISMRMRD::Acquisition >(acq) );
		}
		
		size_t num_acq = all_acquis.items();
		size_t num_iterations = 1000;

		std::vector< PtrAcquisVect > ptr_vectors;

		for(size_t i=0; i<num_iterations; i++)
		{
			std::cout << " NUM IT "<<i<<std::endl;
			PtrAcquisVect temp_sptr_vector;

			for(size_t i_acq=0; i_acq<num_acq; i_acq++)
			{	
				auto sptr_acq = sptr_vec.get_sptr_acq(i_acq);
				temp_sptr_vector.append_sptr_acq(sptr_acq);
			}
			ptr_vectors.push_back( temp_sptr_vector );
		}

		// for(size_t i=0; i<num_iterations; i++)
		// {
		// 	std::cout << "loopindex i " << i <<std::endl;
		// 	AcquisitionsVector temp_dummy_vector;
		// 	temp_dummy_vector.copy_acquisitions_info(all_acquis);
			

		// 	for(size_t i_acq=0; i_acq<all_acquis.number(); i_acq++)
		// 	{	
		// 		ISMRMRD::Acquisition acq;
		// 	    all_acquis.get_acquisition(i_acq, acq);
		// 	    temp_dummy_vector.append_acquisition(acq);
		// 	}

		// 	ISMRMRD::Acquisition acq;	
		// 	for(size_t j=0; j<temp_dummy_vector.items(); j++)
		// 	{
				
		// 	    temp_dummy_vector.get_acquisition(j, acq);
		// 	    std::cout << acq.scan_counter() << std::endl;
		// 	}

		// }

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
		return true;
	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}

