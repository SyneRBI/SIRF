/* ================================================

Author: Johannes Mayer
Date: 2018.04.10
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include "tests_mr_acquisition_model.h"

#include <memory>

#include <ismrmrd/ismrmrd.h>

#include "gadgetron_image_wrap.h"



#include "auxiliary_testing_functions.h"



bool tests_mracqmod::test_fwd_method( void ) 
{
	ISMRMRD::Image< complex_float_t > img = aux_test::get_mock_ismrmrd_image_with_cube();
	CoilDataAsCFImage csm = aux_test::get_mock_coildata_as_cfimage();
	
	ImageWrap img_wrap(MOCK_IMAGE_TYPE, new ISMRMRD::Image< complex_float_t >(img));

	return false;
}



