/* ================================================

Author: Johannes Mayer
Date: 2018.04.03
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tests_auxiliary_testing_functions.h"

#include <ismrmrd/ismrmrd.h>

#include "gadgetron_data_containers.h" 


bool test_aux_test_funs::test_get_mock_acquisition_vector( void )
{

	ISMRMRD::IsmrmrdHeader hdr = aux_test::get_mock_ismrmrd_header();
	AcquisitionsVector acq_vec = aux_test::get_mock_acquisition_vector( hdr );

	return true;

}

