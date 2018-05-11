/* ================================================

Author: Johannes Mayer
Date: 2018.04.10
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include "tests_mr_acquisition_model.h"

#include <memory>
#include <exception>

#include <ismrmrd/ismrmrd.h>

#include "gadgetron_x.h"
#include "gadgetron_image_wrap.h"
#include "gadgetron_data_containers.h"


#include "auxiliary_testing_functions.h"




bool tests_mracqmod::test_fwd_method( void ) 
{
	try
	{
		
		ISMRMRD::Image< complex_float_t > img = aux_test::get_mock_ismrmrd_image_with_cube();
		CoilDataAsCFImage csm = aux_test::get_mock_coildata_as_cfimage();

		ImagesVector img_vec;
		img_vec.append(MOCK_IMAGE_TYPE, new ISMRMRD::Image< complex_float_t> (img));

		ISMRMRD::IsmrmrdHeader hdr = aux_test::get_mock_ismrmrd_header();
		AcquisitionsVector acq_vec = aux_test::get_mock_acquisition_vector( hdr );

		
		MRAcquisitionModel acq_model(std::shared_ptr<AcquisitionsVector> ( new AcquisitionsVector(acq_vec) ), 
							  std::shared_ptr<ImagesVector> (new ImagesVector(img_vec) ));


		ImageWrap img_wrap(MOCK_IMAGE_TYPE, new ISMRMRD::Image< complex_float_t >(img));		
		AcquisitionsVector target_acqs;

		unsigned int offset = 0;
		acq_model.fwd(img_wrap, csm, target_acqs, offset);
		
		unsigned int const num_stored_acquisitions = target_acqs.items();
		unsigned int const num_expected_acquisitions = acq_vec.items();

		std::cout << epiph(num_stored_acquisitions) << std::endl;
		std::cout << epiph(num_expected_acquisitions) << std::endl;


		bool number_acquisitions_match = (num_stored_acquisitions == num_expected_acquisitions);

		target_acqs.write(ACQU_FILE_NAME);

		return number_acquisitions_match;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}


}




