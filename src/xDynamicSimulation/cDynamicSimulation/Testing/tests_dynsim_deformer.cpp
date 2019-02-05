/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include <sstream>

#include "tests_dynsim_deformer.h"

#include "auxiliary_testing_functions.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "sirf/cReg/NiftiImageData3DDeformation.h"



using namespace sirf;

bool DynSimDeformerTester::test_deform_contrast_generator( void )
{
	
try
	{

		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

		MRContrastGenerator mr_cont_gen( segmentation_labels, XML_XCAT_PATH);

		ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
		mr_cont_gen.set_rawdata_header(hdr);



		int const num_simul_cardiac_states = 10;
		MRMotionDynamic  motion_dyn(num_simul_cardiac_states);

		// auto motion_fields = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		// motion_dyn.set_displacement_fields( motion_fields, false );
		
		auto motion_fields = read_cardiac_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		motion_dyn.set_displacement_fields( motion_fields, true );
		
		motion_dyn.prep_displacements_fields();

		mr_cont_gen.map_contrast();
		auto static_state = mr_cont_gen.get_contrast_filled_volumes();

		std::string filename_static = std::string(SHARED_FOLDER_PATH) + "mr_contrast_map_static";
		data_io::write_ISMRMRD_Image_to_Analyze(filename_static, static_state[0]);



		for( size_t i_motion=0; i_motion<num_simul_cardiac_states; i_motion++)
		{
			float const signal = (float)i_motion/(float)num_simul_cardiac_states;

			std::cout << "Getting mvf for signal = " << signal << std::endl;
			NiftiImageData3DDeformation<float> curr_mvf = motion_dyn.get_interpolated_displacement_field( signal );

			std::vector< NiftiImageData3DDeformation<float> > vec_mvfs;
			vec_mvfs.push_back( curr_mvf );

			mr_cont_gen.map_contrast();
			DynamicSimulationDeformer::deform_contrast_generator(mr_cont_gen, vec_mvfs);
			
			auto curr_motion_state = mr_cont_gen.get_contrast_filled_volumes();
			
			std::stringstream filename_stream;
			filename_stream << SHARED_FOLDER_PATH << "mr_contrast_map_state_" << i_motion; 		
			
			data_io::write_ISMRMRD_Image_to_Analyze< complex_float_t > (filename_stream.str(), curr_motion_state[0]);
		}

		return true;
		
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}

bool DynSimDeformerTester::test_deform_pet_contrast_generator( void ) 
{

	try{
		PETContrastGenerator pet_cont_gen = aux_test::get_mock_pet_contrast_generator();

			
		int const num_simul_cardiac_states = 10;
		PETMotionDynamic  motion_dyn(num_simul_cardiac_states);

		// auto motion_fields = read_respiratory_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		// motion_dyn.set_displacement_fields( motion_fields, false );
		
		auto motion_fields = read_cardiac_motionfield_from_h5( H5_XCAT_PHANTOM_PATH );
		motion_dyn.set_displacement_fields( motion_fields, true );
		
		motion_dyn.prep_displacements_fields();

		STIRImageData template_img( PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH ); 
		motion_dyn.align_motion_fields_with_image( template_img );

		pet_cont_gen.map_tissue();
		std::vector< sirf::STIRImageData > static_state = pet_cont_gen.get_contrast_filled_volumes();

		std::string filename_static = std::string(SHARED_FOLDER_PATH) + "/pet_activity_map_static";
		data_io::write_PET_image_to_hv(filename_static, static_state[0]);



		for( size_t i_motion=0; i_motion<num_simul_cardiac_states; i_motion++)
		{
			float const signal = (float)i_motion/(float)num_simul_cardiac_states;

			std::cout << "Getting mvf for signal = " << signal << std::endl;
			NiftiImageData3DDeformation<float> curr_mvf = motion_dyn.get_interpolated_displacement_field( signal );

			std::vector< NiftiImageData3DDeformation<float> > vec_mvfs;
			vec_mvfs.push_back( curr_mvf );

			pet_cont_gen.map_tissue();
			DynamicSimulationDeformer::deform_contrast_generator(pet_cont_gen, vec_mvfs);
			
			std::vector< sirf::STIRImageData > curr_motion_state = pet_cont_gen.get_contrast_filled_volumes();
			
			std::stringstream filename_stream;
			filename_stream << SHARED_FOLDER_PATH << "pet_activity_map_state_" << i_motion; 		
			data_io::write_PET_image_to_hv(filename_stream.str(), curr_motion_state[0]);
		
		}

		return true;
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
	
		typedef NiftiImageData<float> ImageType;
		size_t const num_cycles = 10000;

		

		for(size_t i_cycle=0; i_cycle<num_cycles; i_cycle++)
		{
			ImageType temp_img(DISPLACEMENT_FIELD_PATH);				
			std::cout << "Cycle: " << i_cycle+1 << "/" << num_cycles << std::endl;
			// ImageType another_img(temp_img);				
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