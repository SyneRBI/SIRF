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

		SIRFImageDataDeformation img_dat_deform;

		ISMRMRD::NDArray< unsigned int > segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );
		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		mr_cont_gen.map_contrast();

		// DynamicSimulationDeformer::deform_contrast_generator(mr_cont_gen, img_dat_deform);

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