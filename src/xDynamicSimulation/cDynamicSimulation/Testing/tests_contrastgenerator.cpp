/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tests_contrastgenerator.h"




// contrast generator

bool test_contgen::test_constructor( void )
{

	LabelArray label_arr = aux_test::get_mock_label_array();
	MRContrastGenerator mr_contgen (label_arr, XML_TEST_PATH); 



	return false;
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
