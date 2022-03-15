/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include <sstream>

#include "tests_dynsim_deformer.h"

#include "auxiliary_testing_functions.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "sirf/Reg/NiftyResample.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"

#include "sirf/Gadgetron/mrtest_auxiliary_funs.h"


using namespace sirf;


bool DynSimDeformerTester::test_nifti_data_deformation( void )
{
	std::cout << " --- Running " << __FUNCTION__ << std::endl;
	try
	{
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

		auto sptr_img_to_deform = std::make_shared< NiftiImageData3D<float> >( segmentation_labels );

		// auto motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		auto motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

		for(size_t i=0; i<motion_fields.size(); i++)
		{

			std::stringstream sstream_output;	
			sstream_output << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__ << "_state_" << i; 		

			NiftyResampler<float> resampler;

			resampler.set_interpolation_type_to_cubic_spline();
			resampler.set_reference_image(sptr_img_to_deform);
			resampler.set_floating_image (sptr_img_to_deform);

			auto disp_trafo = std::make_shared<NiftiImageData3DDeformation<float> >( ( motion_fields[i] ).get_as_deformation_field(motion_fields[i]) );
			resampler.add_transformation(disp_trafo);

			resampler.process();

			auto output_img = resampler.get_output_sptr();
			
			output_img->write( sstream_output.str() );
			

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

bool DynSimDeformerTester::test_mr_geometry( void )
{
	std::cout << " --- Running " << __FUNCTION__ << std::endl;

	try
		{
			std::string root_path = "/home/sirfuser/devel/install/share/SIRF-3.1/Simulation/Input/";
			std::string fname_segmentation = root_path + "segmentation.nii";
			std::string fname_ct = root_path + "contrast_template.h5";
			std::string fname_at = root_path + "acquisition_template.h5";

			LabelVolume segmentation = NiftiImageData3D<float>(fname_segmentation);

			segmentation.get_geom_info_sptr()->print_info();

			AcquisitionsVector contrast_rawdata(fname_ct);
			ISMRMRD::MatrixSize matsize = contrast_rawdata.acquisitions_info().get_IsmrmrdHeader().encoding[0].reconSpace.matrixSize;

			int const Nx = matsize.x;
			int const Ny = matsize.y;
			int const Nz = matsize.z;

			GadgetronImagesVector contrast_template(contrast_rawdata);

			for(int i=0; i<contrast_template.number(); ++i)
			{
				auto ptr_img = (CFImage*) contrast_template.sptr_image_wrap(i)->ptr_image();
				ptr_img->resize(Nx, Ny, Nz, 1);
				
				// #pragma omp parallel
				for( size_t nz=0; nz<Nz; nz++)
				for( size_t ny=0; ny<Ny; ny++)
				for( size_t nx=0; nx<Nx; nx++)
				{
					size_t linear_index_access = (nz*Ny + ny)*Nx + nx;
					ptr_img->operator()(nx, ny, nz, 0) = segmentation(nx,ny,nz);						
				}
				
			}
			contrast_template.get_geom_info_sptr()->print_info();
			contrast_template.reorient(*(segmentation.get_geom_info_sptr()));

			std::stringstream name_stream;	
			name_stream << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_contrast" << __FUNCTION__; 		

			sirf::write_imagevector_to_raw(name_stream.str(), contrast_template);

			NiftyResampler<float> resampler;
		    resampler.set_interpolation_type_to_cubic_spline();

			// both floating and reference image must be in the same coordinate system
			resampler.set_floating_image (std::make_shared< GadgetronImagesVector> (contrast_template));
			
			AcquisitionsVector target_rawdata(fname_at);
			GadgetronImagesVector acquisition_template(target_rawdata);
			resampler.set_reference_image(std::make_shared< GadgetronImagesVector> (acquisition_template));

			acquisition_template.get_geom_info_sptr()->print_info();
			contrast_template.get_geom_info_sptr()->print_info();


			resampler.process();

			const std::shared_ptr<const sirf::ImageData>  sptr_deformed_img = resampler.get_output_sptr();

			GadgetronImagesVector img_data(acquisition_template);

			sptr_deformed_img->copy(sptr_deformed_img->begin(),
									img_data.begin(), 
									img_data.end());


			name_stream.str("");
			name_stream << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_resampled" << __FUNCTION__; 		
			sirf::write_imagevector_to_raw(name_stream.str(), img_data);

			return true;
			
		}
		catch( std::runtime_error const &e)
		{
			std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
			std::cout << e.what() << std::endl;
			throw e;
		}
}

bool DynSimDeformerTester::test_deform_mr_contrast_generator( void )
{
	std::cout << " --- Running " << __FUNCTION__ << std::endl;

	try
		{
			LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

			MRContrastGenerator mr_cont_gen (segmentation_labels, XML_TEST_PATH);  	
			sirf::AcquisitionsVector av(ISMRMRD_H5_TEST_PATH);
			mr_cont_gen.set_template_rawdata(av);
			
			auto motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			

			for( size_t i=0; i<motion_fields.size(); i++)
			{
				std::cout << "Getting mvf #" << i << std::endl;
				NiftiImageData3DDeformation<float> curr_mvf = (motion_fields[i]).get_as_deformation_field( motion_fields[i] );

				std::vector< NiftiImageData3DDeformation<float> > vec_mvfs;
				vec_mvfs.push_back( curr_mvf );

				mr_cont_gen.map_contrast();
				DynamicSimulationDeformer dsd;
				dsd.deform_contrast_generator(mr_cont_gen, vec_mvfs);
				
				GadgetronImagesVector curr_motion_state = mr_cont_gen.get_contrast_filled_volumes();

				std::stringstream name_stream;	
				name_stream << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__ << "_state_" << i; 		

				sirf::write_imagevector_to_raw(name_stream.str(), curr_motion_state);

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

	std::cout << " --- Running " << __FUNCTION__ << std::endl;

	try{
		PETContrastGenerator pet_cont_gen = aux_test::get_mock_pet_contrast_generator();

			
		int const num_simul_cardiac_states = 8;
		PETMotionDynamic  motion_dyn(num_simul_cardiac_states);

		auto motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		// auto motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
		
		motion_dyn.set_displacement_fields( motion_fields, true );
		motion_dyn.prep_displacement_fields();

		STIRImageData template_img( PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH ); 

		pet_cont_gen.map_tissue();
		std::vector< sirf::STIRImageData > static_state = pet_cont_gen.get_contrast_filled_volumes();


		for( size_t i_motion=0; i_motion<num_simul_cardiac_states; i_motion++)
		{
			float const signal = (float)i_motion/(float)num_simul_cardiac_states;

			std::cout << "Getting mvf for signal = " << signal << std::endl;
			NiftiImageData3DDeformation<float> curr_mvf = motion_dyn.get_interpolated_deformation_field( signal );

			std::vector< NiftiImageData3DDeformation<float> > vec_mvfs;
			vec_mvfs.push_back( curr_mvf );

			pet_cont_gen.map_tissue();
			DynamicSimulationDeformer dsd;
			dsd.deform_contrast_generator(pet_cont_gen, vec_mvfs);
			
			std::vector< sirf::STIRImageData > curr_motion_state = pet_cont_gen.get_contrast_filled_volumes();
			
			std::stringstream filename_stream;	
			filename_stream << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__ << "_state_" << i_motion; 		


			sirf::NiftiImageData3D<float> curr_motion_state_nii( curr_motion_state[0] );				
			curr_motion_state_nii.write(filename_stream.str() );
			// data_io::write_PET_image_to_hv(filename_stream.str(), curr_motion_state[0]);		
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

bool DynSimDeformerTester::test_motion_of_MotionDynamics()
{
	std::cout << " --- Running " << __FUNCTION__ << std::endl;

	try
	{
		bool test_succesful = true;

		//	---------------------------------------
		// Prep
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

		MRContrastGenerator mr_cont_gen (segmentation_labels, XML_TEST_PATH);  	
		sirf::AcquisitionsVector av(ISMRMRD_H5_TEST_PATH);
		mr_cont_gen.set_template_rawdata(av);

		//	---------------------------------------
		// Set motion infos
		int const num_simul_motion_dyn = 8;
		MRMotionDynamic motion_dyn( num_simul_motion_dyn );
		
		bool card_motion = false;

		if( card_motion )
		{
			auto motion_fields = read_cardiac_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			motion_dyn.set_displacement_fields( motion_fields, card_motion );
		}
		else
		{
			auto motion_fields = read_respiratory_motionfields_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
			motion_dyn.set_displacement_fields( motion_fields, card_motion );
		}

 		
		std::vector< SignalBin > motion_bins = motion_dyn.get_bins();

		for( size_t i=0; i<motion_bins.size(); i++)
		{

			float const motion_state = std::get<1>(motion_bins[i]);

			std::cout << "Getting mvf for motion state " << motion_state << std::endl;
			
			sirf::NiftiImageData3DDeformation<float> curr_mvf = motion_dyn.get_interpolated_deformation_field( motion_state );
			
			std::vector< NiftiImageData3DDeformation<float> > vec_mvfs;
			vec_mvfs.push_back( curr_mvf );

			mr_cont_gen.map_contrast();

			DynamicSimulationDeformer dsd;
			dsd.deform_contrast_generator(mr_cont_gen, vec_mvfs);
			
			GadgetronImagesVector curr_motion_state = mr_cont_gen.get_contrast_filled_volumes();
			std::stringstream name_stream;
			name_stream << SHARED_FOLDER_PATH << TESTDATA_OUT_PREFIX << "output_" << __FUNCTION__ << "_dyn_" << i; 		
			sirf::write_imagevector_to_raw(name_stream.str(), curr_motion_state);

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