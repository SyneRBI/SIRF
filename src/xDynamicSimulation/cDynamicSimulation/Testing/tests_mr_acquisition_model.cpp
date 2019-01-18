/* ================================================

Author: Johannes Mayer
Date: 2018.04.10
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include "tests_mr_acquisition_model.h"

#include <memory>
#include <exception>

#include <ismrmrd/ismrmrd.h>

#include "sirf/cGadgetron/gadgetron_x.h"
#include "sirf/cGadgetron/gadgetron_image_wrap.h"
#include "sirf/cGadgetron/gadgetron_data_containers.h"


#include "auxiliary_testing_functions.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

using namespace sirf;


bool tests_mracqmod::test_fwd_method( void ) 
{
	try
	{
		
		ISMRMRD::Image< complex_float_t > img = aux_test::get_mock_ismrmrd_image_with_cube();

		std::vector<float> img_data;
		img_data.resize(img.getNumberOfDataElements(), 0);
		for (int i_pixel=0; i_pixel<img.getNumberOfDataElements(); i_pixel++)
		{
			img_data[i_pixel] = std::real( *(img.begin() + i_pixel));			
		}

		data_io::write_raw<float> ( std::string(SHARED_FOLDER_PATH) + "test_fwd_method_imgInput", &img_data[0], img_data.size());

		CoilDataAsCFImage csm = aux_test::get_mock_coildata_as_cfimage();

        sirf::GadgetronImagesVector img_vec;
		img_vec.append(MOCK_DATA_TYPE, new ISMRMRD::Image< complex_float_t> (img));

		// ISMRMRD::IsmrmrdHeader hdr = aux_test::get_mock_ismrmrd_header();
		// AcquisitionsVector source_acqs = aux_test::get_mock_acquisition_vector( hdr );
		
		ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
		AcquisitionsVector source_acqs = mr_io::read_ismrmrd_acquisitions(ISMRMRD_H5_TEST_PATH);
		
		
		MRAcquisitionModel acq_model(std::shared_ptr<AcquisitionsVector> ( new AcquisitionsVector(source_acqs) ), 
							  std::shared_ptr<GadgetronImagesVector> (new GadgetronImagesVector(img_vec) ));

		ImageWrap img_wrap(MOCK_DATA_TYPE, new ISMRMRD::Image< complex_float_t >(img));		

		AcquisitionsVector target_acqs;
		target_acqs.copy_acquisitions_info( source_acqs );

        complex_float_t *img_wrap_data = new complex_float_t[img.getNumberOfDataElements()];
		
        img_wrap.get_complex_data(img_wrap_data);
        img_wrap_data->real();

        float *re_img_wrap_data = new float[img.getNumberOfDataElements()];
        for (int i=0; i<img.getNumberOfDataElements(); ++i)
            re_img_wrap_data[i] = img_wrap_data[i].real();
		
		data_io::write_raw<float> ( std::string(SHARED_FOLDER_PATH) + "test_fwd_method_real_imgWrapData", re_img_wrap_data, img.getNumberOfDataElements());

        delete [] img_wrap_data;
        img_wrap_data = nullptr;
        delete [] re_img_wrap_data;
        re_img_wrap_data = nullptr;
	
		unsigned int offset = 0;

		acq_model.fwd(img_wrap, csm, target_acqs, offset);
		
		unsigned int const num_stored_acquisitions = target_acqs.items();
		unsigned int const num_expected_acquisitions = source_acqs.items();

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




