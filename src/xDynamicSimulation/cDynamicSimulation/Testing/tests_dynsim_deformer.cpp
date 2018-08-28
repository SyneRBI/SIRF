/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "auxiliary_testing_functions.h"

#include "tests_dynsim_deformer.h"


bool DynSimDeformerTester::test_deform_contrast_generator( void )
{
	
try
	{
		bool test_succesful = true;

		ISMRMRD::NDArray< unsigned int > segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );

		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
		mr_cont_gen.set_rawdata_header(hdr);

		mr_cont_gen.map_contrast();

		SIRFImageDataDeformation img_dat_deform( DISPLACEMENT_FIELD_PATH );


		auto sptr_mvf_as_nifti = img_dat_deform.get_image_as_nifti();
		nifti_image mvf_as_nifti = *sptr_mvf_as_nifti;


		std:: cout << epiph(mvf_as_nifti.ndim) << std::endl;
		std:: cout << epiph(mvf_as_nifti.nx) << std::endl;
		std:: cout << epiph(mvf_as_nifti.ny) << std::endl;
		std:: cout << epiph(mvf_as_nifti.nz) << std::endl;
		std:: cout << epiph(mvf_as_nifti.nt) << std::endl;
		std:: cout << epiph(mvf_as_nifti.nu) << std::endl;



		DynamicSimulationDeformer::deform_contrast_generator(mr_cont_gen, img_dat_deform);

		test_succesful = false;

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}




}