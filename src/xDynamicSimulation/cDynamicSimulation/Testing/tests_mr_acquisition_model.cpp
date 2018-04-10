/* ================================================

Author: Johannes Mayer
Date: 2018.04.10
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include "tests_mr_acquisition_model.h"

#include <memory>

#include <ismrmrd/ismrmrd.h>

#include "gadgetron_x.h"
#include "gadgetron_image_wrap.h"



#include "auxiliary_testing_functions.h"



void tests_mracqmod::test_get_serialized_ismrmrd_header( void )
{

	std::string serialized_hdr = aux_test::get_serialized_mock_ismrmrd_header();

	std::cout << serialized_hdr << std::endl;

}








bool tests_mracqmod::test_fwd_method( void ) 
{
	
	ISMRMRD::Image< complex_float_t > img = aux_test::get_mock_ismrmrd_image_with_cube();
	CoilDataAsCFImage csm = aux_test::get_mock_coildata_as_cfimage();
	

	ImagesVector img_vec;
	img_vec.append(MOCK_IMAGE_TYPE, new ISMRMRD::Image< complex_float_t> (img));
	
	ISMRMRD::Acquisition acq = aux_test::get_mock_ismrmrd_acquisition();

	AcquisitionsVector acq_vec;
	acq_vec.append_acquisition(acq);


	MRAcquisitionModel ma(std::shared_ptr<AcquisitionsVector> ( new AcquisitionsVector(acq_vec) ), 
						  std::shared_ptr<ImagesVector> (new ImagesVector(img_vec) ));
		

	ImageWrap img_wrap(MOCK_IMAGE_TYPE, new ISMRMRD::Image< complex_float_t >(img));		
	AcquisitionsVector target_acqs;

	unsigned int offset = 0;
	ma.fwd(img_wrap, csm, target_acqs, offset);


	return false;
}




