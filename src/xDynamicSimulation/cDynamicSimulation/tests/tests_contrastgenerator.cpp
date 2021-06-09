/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tests_contrastgenerator.h"


#include <string>
#include <sstream>
#include <stdio.h>
#include <iostream>

#include <ismrmrd/xml.h>

#include "auxiliary_testing_functions.h"

#include "sirf/cDynamicSimulation/tissuelabelmapper.h"
#include "sirf/cDynamicSimulation/tissueparameters.h"
#include "sirf/cDynamicSimulation/contrastgenerator.h"
#include "sirf/cDynamicSimulation/phantom_input.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "mrtest_auxiliary_funs.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"


using namespace std;
using namespace ISMRMRD;
using namespace stir;
using namespace sirf;

// contrast generator


bool test_contgen::test_get_tissue_parameter( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{	
		MRContrastGenerator cont_gen = aux_test::get_mock_mr_contrast_generator();

		LabelType label = 52;
		TissueParameter tiss_par = cont_gen.get_petmr_tissue_parameter( label );		

		std::cout << "Label " << label << " was associated with the parameter " << tiss_par.name_ <<std::endl;

		return true;
	}
	catch( const std::exception& e)
	{	
		std::cout << "Exception caught at highest level in main" << std::endl;
		std::cout<< e.what() << '\n';			
		return false;
	}

}


bool test_contgen::test_mr_constructor( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	LabelVolume segmentation_labels = aux_test::get_mock_label_volume();
	MRContrastGenerator mr_contgen (segmentation_labels, XML_TEST_PATH); 

	return true;
}


bool test_contgen::test_mr_set_rawdata_header( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{	
		LabelVolume segmentation_labels = aux_test::get_mock_label_volume();
		MRContrastGenerator mr_contgen (segmentation_labels, XML_TEST_PATH);  	

		ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
		mr_contgen.set_rawdata_header(hdr);

		return  true;
	}
	catch(...)
	{
		std::cout << "An unknown exception was caught" << std::endl;
		return false;
	}
}



bool test_contgen::test_mr_map_contrast_dim_check( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	LabelVolume segmentation_labels = aux_test::get_mock_label_volume();
	MRContrastGenerator mr_contgen (segmentation_labels, XML_TEST_PATH);  	
	ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);

	mr_contgen.set_rawdata_header(hdr);
	mr_contgen.map_contrast();

	GadgetronImagesVector& contrasts = mr_contgen.get_contrast_filled_volumes();	


	size_t const num_echoes = ((hdr.sequenceParameters.get()).TE.get()).size();;
	size_t input_dims[8] = {3,MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE,num_echoes,1,1,1};

	std::vector< int > contrast_dims;

	sirf::Dimensions dims = contrasts.dimensions();

	contrast_dims.push_back( 3 );
	contrast_dims.push_back( dims["x"] );
	contrast_dims.push_back( dims["y"] );
	contrast_dims.push_back( dims["z"] );
	contrast_dims.push_back( dims["n"] );

	bool dims_are_correct = true; 
	for( int i=0; i< 4; i++)
		dims_are_correct *= (contrast_dims[i] == input_dims[i]);
		
	return dims_are_correct;
}

void test_contgen::test_mr_map_contrast_application_to_xcat( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;
	
	std::cout << "Reading segmentation ... " <<std::endl;
	LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
	std::cout << "... finished. " <<std::endl;
	
	MRContrastGenerator mr_contgen( segmentation_labels, XML_XCAT_PATH);
	ISMRMRD::IsmrmrdHeader hdr =  mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
	mr_contgen.set_rawdata_header(hdr);

	
	mr_contgen.map_contrast();
	
	GadgetronImagesVector& mr_contrasts = mr_contgen.get_contrast_filled_volumes();	
	
	std::stringstream name_stream;
	name_stream << SHARED_FOLDER_PATH << "testoutput_mr_cont_gent_contrast_";
	sirf::write_imagevector_to_raw(name_stream.str(), mr_contrasts);
}

void test_contgen::test_get_signal_for_tissuelabel_in_xcat()
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );
	MRContrastGenerator mr_contgen( segmentation_labels, XML_XCAT_PATH);
	IsmrmrdHeader hdr =  mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);

	mr_contgen.set_rawdata_header(hdr);


	size_t const test_label = 3;

	auto signal = mr_contgen.get_signal_for_tissuelabel( test_label );

	std::cout <<"Signal in label " << test_label << " amounts to: " << signal << std::endl;

}

void test_contgen::test_replace_petmr_tissue_parameters_in_xcat()
{

	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

	MRContrastGenerator mr_contgen( segmentation_labels, XML_XCAT_PATH);
	IsmrmrdHeader hdr =  mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
	
	mr_contgen.set_rawdata_header(hdr);
	mr_contgen.map_contrast();

	GadgetronImagesVector& mr_contrasts = mr_contgen.get_contrast_filled_volumes();
	std::string output_name = std::string(SHARED_FOLDER_PATH) + "test_contrast_gen_pre_contrast";
	sirf::write_imagevector_to_raw(output_name, mr_contrasts);

	// now replace one label and see what happens in the image
	LabelType label_to_replace = 1;
	auto tissue_param_pair = aux_test::get_mock_contrast_signal_extremes();

	mr_contgen.replace_petmr_tissue_parameters(label_to_replace, tissue_param_pair.second);
	mr_contgen.map_contrast();

	mr_contrasts = mr_contgen.get_contrast_filled_volumes();
	output_name = std::string(SHARED_FOLDER_PATH) + "test_contrast_gen_post_contrast";
	sirf::write_imagevector_to_raw(output_name, mr_contrasts);
}



bool test_contgen::test_map_flash_contrast( void )
{

	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	TissueParameter tiss_par = aux_test::get_mock_tissue_parameter();
	auto ptr_to_mock_tiss = std::make_shared<TissueParameter>(tiss_par);

	ISMRMRD::IsmrmrdHeader hdr = aux_test::get_mock_ismrmrd_header();

	std::vector <complex_float_t> flash_contrast = map_flash_contrast(ptr_to_mock_tiss, hdr);


	float const t1 = 1;
	float const t2 = 2;
	float const dens = 100;
	float const angle = M_PI/2;
	float const cs = 1;
	float const field_strength_t = 1.0;

	float const TR = 2;
	float const TE = 1;


	complex_float_t IMAG_UNIT(0,1);

	complex_float_t input_contrast_echo1 = exp( IMAG_UNIT * (float)42.58/1000.f * TE * cs * field_strength_t)*dens * (float)sin(angle) * 
	(float)(1-exp(-TR/t1)) / (float)(1- exp(-TR/t1)*cos(angle)) * (float)exp(-TE/t2);	
	complex_float_t mock_contrast = flash_contrast[0];

	float const epsilon = 0.00001;

	bool equal_contrast = (input_contrast_echo1.real() - mock_contrast.real() ) < epsilon ;
	equal_contrast *= (input_contrast_echo1.imag() - mock_contrast.imag() < epsilon );

	return equal_contrast;

}

bool test_contgen::test_pet_constructor( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{
		LabelVolume segmentation_labels = aux_test::get_mock_label_volume();
		PETContrastGenerator pet_contgen (segmentation_labels, XML_XCAT_PATH); 

		return true;
	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}

bool test_contgen::test_pet_map_contrast( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

		PETContrastGenerator pet_contgen (segmentation_labels, XML_XCAT_PATH); 
		pet_contgen.set_template_image_from_file( std::string(PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH) );						

		pet_contgen.map_contrast();



		return true;
	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}


bool test_contgen::test_pet_map_attenuation( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

		PETContrastGenerator pet_contgen( segmentation_labels, XML_XCAT_PATH ); 
		pet_contgen.set_template_image_from_file( PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH );						

		pet_contgen.map_attenuation();

		pet_contgen.map_contrast();

		return true;
	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
}

bool test_contgen::test_set_template_image_from_file( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

		PETContrastGenerator pet_contgen( segmentation_labels, XML_XCAT_PATH ); 
		pet_contgen.set_template_image_from_file( PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH );						


		auto voxel_sizes = pet_contgen.get_voxel_sizes();
		auto dims = pet_contgen.get_dimensions();

		for(int i=0;i<3;i++)
		{
			std::cout<< dims[i] << std::endl;	
			std::cout<< voxel_sizes[i] << std::endl;
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

bool test_contgen::test_resample_to_template_image( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{
		LabelVolume segmentation_labels = read_segmentation_to_nifti_from_h5( H5_XCAT_PHANTOM_PATH );

		PETContrastGenerator pet_contgen( segmentation_labels, XML_XCAT_PATH ); 
		pet_contgen.set_template_image_from_file( PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH );						

		pet_contgen.map_contrast();

		auto resampled_images = pet_contgen.get_contrast_filled_volumes( true );

		std::cout << "There are " << resampled_images.size() << " resampled PET images." <<std::endl;

		for(int i=0; i<resampled_images.size(); ++i)
		{
			stringstream outname;
			outname << SHARED_FOLDER_PATH << "resampled_PET_image_" << i;

			sirf::NiftiImageData3D<float> resampled_nifti( resampled_images[i] );

			resampled_nifti.write( outname.str() );
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

void test_contgen::test_pet_map_contrast_application_to_xcat( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{
		PETContrastGenerator pet_contgen = aux_test::get_mock_pet_contrast_generator( );

		pet_contgen.map_contrast();
		auto volume_container = pet_contgen.get_contrast_filled_volumes();

		STIRImageData contrast_volume = volume_container[0];
	
		auto dims = pet_contgen.get_dimensions();

		int Nx = dims[0];
		int Ny = dims[1];
		int Nz = dims[2];
		std::cout << " dims " << Nx << "," << Ny << ","<< Nz;
		std::stringstream outname_contrast; 
		outname_contrast << std::string(SHARED_FOLDER_PATH) << "xcat_pet_contrast";


		std::vector<float> data_output(Nx*Ny*Nz, 0.f);
		// for( size_t i=0; i<data_output.size(); i++)
			// std::cout << data_output[i] <<std::endl;

		contrast_volume.get_data( &data_output[0] );

		// contrast_volume.write( outname_contrast.str() );
		// for( size_t i=0; i<data_output.size(); i++)
			// std::cout << data_output[i] <<std::endl;

		data_io::write_raw< float >( outname_contrast.str() , &data_output[0], Nx*Ny*Nz);

	}
	catch( std::runtime_error const &e)
	{	
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
	
}


// tissue label mapper 

bool test_tlm::test_get_filepath_tissue_parameter_xml()
{

	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	LabelVolume labels_list = aux_test::get_mock_label_volume();
	TissueLabelMapper tlm(labels_list, XML_TEST_PATH);

	std::string xml_path = tlm.get_filepath_tissue_parameter_xml();

	if (xml_path.compare(XML_TEST_PATH) == 0)
		return true;
}

bool test_tlm::test_get_labels_array()
{
	LabelVolume labels_list = aux_test::get_mock_label_volume();
	TissueLabelMapper tlm(labels_list, XML_TEST_PATH);
	
	LabelVolume reference_list = tlm.get_segmentation_labels();

	return labels_list == reference_list;
}

bool test_tlm::test_get_segmentation_dimensions( void )
{

	LabelVolume labels_list = aux_test::get_mock_label_volume();
	TissueLabelMapper tlm( labels_list, XML_TEST_PATH);

	const int* data_dims = tlm.get_segmentation_dimensions();

	int input_dims[8] = {3,MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE,1,1,1,1};

	bool dims_are_correct = true;

	for( int i=0; i<8; i++)
		dims_are_correct *= (data_dims[i] == input_dims[i]);		
		
	return dims_are_correct;
}


bool test_tlm::test_assign_tissue_parameters_label_found( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	TissueParameterList tiss_list = aux_test::get_mock_tissue_param_list();
	LabelVolume labels_list = aux_test::get_mock_label_volume();

	TissueVector tissue_volume = assign_tissue_parameters_to_labels( tiss_list, labels_list);

	size_t num_elements_tissue_pointers = tissue_volume.size();

	bool all_labels_correct = true;

	for( int i=0; i<num_elements_tissue_pointers; i++)
	{
		
		std::shared_ptr<TissueParameter> current_tissue_param = tissue_volume[i];
		unsigned int associated_label = current_tissue_param->label_;
		
		all_labels_correct *= (labels_list(i) == associated_label);
		
	}

	return all_labels_correct;
}

bool test_tlm::test_assign_tissue_parameters_label_not_found( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	TissueParameterList tiss_list = aux_test::get_mock_tissue_param_list();
	LabelVolume labels_list = aux_test::get_mock_label_volume();
	labels_list(0) = 23;
	try
	{
		TissueVector tissue_volume = assign_tissue_parameters_to_labels( tiss_list, labels_list);
	}
	catch( std::runtime_error const &e)
	{	

		std::cout << "Test output: " << e.what() << std::endl;
		return true;
	}
}

bool test_tlm::test_map_labels_to_tissue_from_xml( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	LabelVolume lab_arr = aux_test::get_mock_label_volume();

	TissueLabelMapper tlm(lab_arr, XML_TEST_PATH);

	tlm.map_labels_to_tissue_from_xml();

	TissueVector tiss_vec = tlm.get_segmentation_tissues();

	bool all_labels_correct = true;
	for (int i = 0; i<tiss_vec.size(); i++)
	{
		unsigned int const tissue_label = tiss_vec[i]->label_;
		
		all_labels_correct *= (tissue_label == lab_arr(i));
	}

	return all_labels_correct;
}


bool test_tlm::test_replace_petmr_tissue_parameters( void )
{
	std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{
		
		LabelVolume lab_arr = aux_test::get_mock_label_volume();
		TissueLabelMapper tlm(lab_arr, XML_TEST_PATH);

		tlm.map_labels_to_tissue_from_xml();

		TissueParameterList tpl_before_replacement = tlm.get_tissue_parameter_list();


		LabelType label_to_replace = 2;

		TissueParameter tiss_par_to_substitute;
		tiss_par_to_substitute.name_ = "lala";
		tiss_par_to_substitute.label_ = 19;
		tiss_par_to_substitute.mr_tissue_.t1_miliseconds_ = 1.01;
		tiss_par_to_substitute.mr_tissue_.t2_miliseconds_ = 1.02;
		tiss_par_to_substitute.mr_tissue_.cs_ppm_ = 1.03;
		tiss_par_to_substitute.mr_tissue_.spin_density_percentH2O_ = 1.04;

		tiss_par_to_substitute.pet_tissue_.attenuation_1_by_cm_= 1.05;
		tiss_par_to_substitute.pet_tissue_.activity_kBq_ml_= 1.06;

		tlm.replace_petmr_tissue_parameters(label_to_replace, tiss_par_to_substitute);

		TissueParameterList tpl_after_replacement = tlm.get_tissue_parameter_list();

		for(size_t i=0; i<tpl_after_replacement.size(); i++)
		{

			TissueParameter curr_tiss = tpl_after_replacement[i];
			TissueParameter ancient_tiss = tpl_before_replacement[i];

			if(curr_tiss.label_ == label_to_replace )	
			{
				std::cout << epiph(curr_tiss.mr_tissue_.t1_miliseconds_) << std::endl;	
				std::cout << epiph(ancient_tiss.mr_tissue_.t1_miliseconds_) << std::endl;
			}
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