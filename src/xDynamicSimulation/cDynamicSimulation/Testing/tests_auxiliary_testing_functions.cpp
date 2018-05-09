/* ================================================

Author: Johannes Mayer
Date: 2018.04.03
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tests_auxiliary_testing_functions.h"

#include <ismrmrd/ismrmrd.h>

#include "gadgetron_data_containers.h" 


bool test_aux_test_funs::test_get_serialized_ismrmrd_header( void )
	{
		try
		{
			std::string serialized_hdr = aux_test::get_serialized_mock_ismrmrd_header();
			std::cout << serialized_hdr << std::endl;
		}
		catch( std::runtime_error const &e)
		{
			std::cout << "Exception caught in " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
		}
		return true;
	}



bool test_aux_test_funs::test_get_mock_acquisition_vector( void )
{

	ISMRMRD::IsmrmrdHeader hdr = aux_test::get_mock_ismrmrd_header();
	AcquisitionsVector acq_vec = aux_test::get_mock_acquisition_vector( hdr );

	unsigned int num_acquis = acq_vec.number();
	std::cout<< epiph(num_acquis) << std::endl;


	ISMRMRD::Acquisition acq;
	int const check_aqu_num[3] = {0, 10, 100};

	for( int i=0; i<3; i++)
	{	
		std::cout << epiph( check_aqu_num[i] ) << std::endl;
		acq_vec.get_acquisition(check_aqu_num[i], acq);


		uint16_t const available_channels = acq.available_channels();
		std::cout << epiph( available_channels ) << std::endl;
		std::cout << epiph( acq.getHead().idx.kspace_encode_step_1 ) << std::endl;
		std::cout << epiph( acq.getHead().idx.kspace_encode_step_2 ) << std::endl;
	}

	std::cout << "\n \n" << std::endl;
	

	for( size_t i=0; i<num_acquis; i++ )
	{
		acq_vec.get_acquisition(i, acq);

		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
		{
			std::cout << epiph( acq.getHead().idx.kspace_encode_step_1 ) << std::endl;
			std::cout << epiph( acq.getHead().idx.kspace_encode_step_2 ) << std::endl;
		}
	}

	return true;

}

