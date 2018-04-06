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

#include "tissuelabelmapper.h"
#include "tissueparameters.h"
#include "contrastgenerator.h"
#include "phantom_input.h"
#include "../auxiliary_input_output.h"


using ISMRMRD::ISMRMRD_NDARRAY_MAXDIM;

// contrast generator

bool test_contgen::test_mr_constructor( void )
{

	LabelArray label_arr = aux_test::get_mock_label_array();
	MRContrastGenerator mr_contgen (label_arr, XML_TEST_PATH); 

	return true;
}


bool test_contgen::test_mr_set_rawdata_header( void )
{
	try
	{	
		LabelArray label_arr = aux_test::get_mock_label_array();
		MRContrastGenerator mr_contgen (label_arr, XML_TEST_PATH);  	

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

 	//using ISMRMRD::ISMRMRD_NDARRAY_MAXDIM;

	LabelArray label_arr = aux_test::get_mock_label_array();
	MRContrastGenerator mr_contgen (label_arr, XML_TEST_PATH);  	

	ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
	mr_contgen.set_rawdata_header(hdr);


	mr_contgen.map_contrast();

	ISMRMRD::NDArray< complex_float_t >	contrast = mr_contgen.get_contrast_filled_volume();	

	int const num_echoes = 3;

	size_t input_dims[ISMRMRD_NDARRAY_MAXDIM] = {2,2,2,num_echoes,1,1,1};
	const size_t* contrast_dims = contrast.getDims();

	bool dims_are_correct = true; 
	for( int i=0; i<ISMRMRD_NDARRAY_MAXDIM; i++)
		dims_are_correct *= (contrast_dims[i] == input_dims[i]);

	return dims_are_correct;
}


void test_contgen::test_mr_map_contrast_application_to_xcat( void )
{
	ISMRMRD::NDArray< unsigned int > segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );

	std::string name_output_segmentation = "/media/sf_SharedFiles/tissue_seg_xcat_test_192x192x192";
	data_io::write_raw<unsigned int>(name_output_segmentation, segmentation_labels.begin(), segmentation_labels.getNumberOfElements());
	

	MRContrastGenerator mr_contgen( segmentation_labels, XML_XCAT_PATH);
	ISMRMRD::IsmrmrdHeader hdr =  mr_io::read_ismrmrd_header(ISMRMRD_H5_TEST_PATH);
	mr_contgen.set_rawdata_header(hdr);

	mr_contgen.map_contrast();

	ISMRMRD::NDArray< complex_float_t >	mr_contrast = mr_contgen.get_contrast_filled_volume();	

	size_t num_elements = mr_contrast.getNumberOfElements();
	std::cout << epiph(num_elements) << std::endl;

	// check data sizes
	const size_t* data_dimension = mr_contrast.getDims();
	std::vector < size_t > dims(data_dimension, data_dimension+ISMRMRD_NDARRAY_MAXDIM);

	ISMRMRD::NDArray< float > mr_contrast_abs, mr_contrast_arg; 
	mr_contrast_abs.resize( dims );
	mr_contrast_arg.resize( dims );

	for( size_t i=0; i<num_elements; i++ )
	{
		*(mr_contrast_abs.begin() + i) = std::abs( *(mr_contrast.begin() + i) );
		*(mr_contrast_arg.begin() + i) = std::arg( *(mr_contrast.begin() + i) );

	}
			
	std::string name_output_contrast  = "/media/sf_SharedFiles/flash_contrast_xcat_test_";

	data_io::write_raw<float>(name_output_contrast + "abs_192x192x192" , mr_contrast_abs.begin(), mr_contrast_abs.getNumberOfElements());
	data_io::write_raw<float>(name_output_contrast + "arg_192x192x192" , mr_contrast_arg.begin(), mr_contrast_arg.getNumberOfElements());
	
}


bool test_contgen::test_map_flash_contrast( void )
{

	TissueParameter tiss_par = aux_test::get_mock_tissue_parameter();
	TissueParameter* const ptr_to_mock_tiss = &tiss_par;

	ISMRMRD::IsmrmrdHeader hdr = aux_test::get_mock_ismrmrd_header();
	ISMRMRD::IsmrmrdHeader* ptr_to_mock_hdr = &hdr;

	std::vector <complex_float_t> flash_contrast = map_flash_contrast(ptr_to_mock_tiss, ptr_to_mock_hdr);


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

	float const epsilon = 0.0000001;

	bool equal_contrast = (input_contrast_echo1.real() - mock_contrast.real() < epsilon );
	equal_contrast *= (input_contrast_echo1.imag() - mock_contrast.imag() < epsilon );

	return equal_contrast;

}




// tissue label mapper 

bool test_tlm::test_get_filepath_tissue_parameter_xml()
{

	LabelArray labels_list = aux_test::get_mock_label_array();
	TissueLabelMapper tlm(labels_list, XML_TEST_PATH);

	std::string xml_path = tlm.get_filepath_tissue_parameter_xml();

	if (xml_path.compare(XML_TEST_PATH) == 0)
		return true;
}

bool test_tlm::test_get_labels_array()
{
	LabelArray labels_list = aux_test::get_mock_label_array();
	TissueLabelMapper tlm(labels_list, XML_TEST_PATH);
	
	LabelArray reference_list = tlm.get_segmentation_labels();

	bool set_and_get_are_the_same = aux_test::equal_array_content<unsigned int> (labels_list, reference_list);

	return set_and_get_are_the_same;
}

bool test_tlm::test_get_segmentation_dimensions( void )
{

	LabelArray labels_list = aux_test::get_mock_label_array();
	TissueLabelMapper tlm( labels_list, XML_TEST_PATH);

	const size_t* data_dims = tlm.get_segmentation_dimensions();

	size_t input_dims[ISMRMRD_NDARRAY_MAXDIM] = {2,2,2,0,0,0,0};

	bool dims_are_correct = true;

	for( int i=0; i<ISMRMRD_NDARRAY_MAXDIM; i++)
		dims_are_correct *= (data_dims[i] == input_dims[i]);		
	
	return dims_are_correct;
}


bool test_tlm::test_assign_tissue_parameters_label_found( void )
{

	TissueParameterList tiss_list = aux_test::get_mock_tissue_param_list();
	LabelArray labels_list = aux_test::get_mock_label_array();

	TissueVector tissue_volume = assign_tissue_parameters_to_labels( tiss_list, labels_list);

	size_t num_elements_tissue_pointers = tissue_volume.size();

	bool all_labels_correct = true;

	for( int i=0; i<num_elements_tissue_pointers; i++)
	{
		
		TissueParameter* current_tissue_param = tissue_volume[i];
		unsigned int associated_label = current_tissue_param->label_;
		
		all_labels_correct *= (labels_list(i) == associated_label);
		
	}

	return all_labels_correct;
}

bool test_tlm::test_assign_tissue_parameters_label_not_found( void )
{

	TissueParameterList tiss_list = aux_test::get_mock_tissue_param_list();
	LabelArray labels_list = aux_test::get_mock_label_array();
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
	LabelArray lab_arr = aux_test::get_mock_label_array();

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
