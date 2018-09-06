/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "auxiliary_testing_functions.h"

#include "SIRFImageDataDeformation.h"

#include "tests_dynsim_deformer.h"


using namespace sirf;

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

		std::vector< SIRFImageDataDeformation > all_img_deforms;
		all_img_deforms.push_back(img_dat_deform);

		auto sptr_mvf_as_nifti = img_dat_deform.get_image_as_nifti();
		nifti_image mvf_as_nifti = *sptr_mvf_as_nifti;

		DynamicSimulationDeformer::deform_contrast_generator(mr_cont_gen, all_img_deforms);

		auto cont_filled_vols = mr_cont_gen.get_contrast_filled_volumes();

		for(int i=0; i<cont_filled_vols.size(); i++)
		{
			std::stringstream output_name;
			output_name << FILENAME_MR_DEFORM_TEST << "_" << i;
			data_io::write_ISMRMRD_Image_to_Analyze< complex_float_t > (output_name.str(), cont_filled_vols[i]);
		}

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}


bool DynSimDeformerTester::test_SIRFImageDataDeformation_memory_behavior()
{

	

	try
	{
		bool test_succesful = true;
	
		typedef SIRFImageData ImageType; 
		size_t const num_cycles = 10000;

		for(size_t i_cycle=0; i_cycle<num_cycles; i_cycle++)
		{
			std::cout << "Cycle: " << i_cycle+1 << "/" << num_cycles << std::endl;
			ImageType temp_img(DISPLACEMENT_FIELD_PATH);				
		}

		return test_succesful;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}

}