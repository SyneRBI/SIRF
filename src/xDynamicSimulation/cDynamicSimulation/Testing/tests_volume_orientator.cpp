/* ================================================

Author: Johannes Mayer
Date: 2018.09.18
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tests_volume_orientator.h"


#include "auxiliary_testing_functions.h"
#include "auxiliary_input_output.h"



#define OUTPUTNAME_INPUT_IMG SHARED_FOLDER_PATH "volume_orientation_input"
#define OUTPUTNAME_REORIENTED_IMG SHARED_FOLDER_PATH "volume_orientation_output"


using namespace sirf;

bool aVolumeOrientatorTester::test_reorient_image(  )
{

	try
	{
		bool test_succesful = true;
		sirf::ReadoutDir ro_dir = ro_dir_z;
		aVolumeOrientator vol_or(ro_dir);

		auto img = aux_test::get_mock_ismrmrd_image_with_gradients();
		data_io::write_ISMRMRD_Image_to_Analyze(OUTPUTNAME_INPUT_IMG, img);


		img = vol_or.reorient_image< float > ( img );
		data_io::write_ISMRMRD_Image_to_Analyze(OUTPUTNAME_REORIENTED_IMG, img);

		return test_succesful;

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}

