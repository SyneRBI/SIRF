/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once




namespace test_contgen
{

bool test_mr_constructor( void );
bool test_mr_set_rawdata_header( void );

bool test_map_flash_contrast( void );

bool test_mr_map_contrast_dim_check( void );

void test_match_output_dims_to_headerinfo( void );

void test_mr_map_contrast_application_to_xcat( void );


bool test_pet_constructor( void );
bool test_pet_map_contrast( void );
bool test_pet_map_attenuation( void );
bool test_set_template_image_from_file( void );

void test_pet_map_contrast_application_to_xcat( void );
void test_replace_petmr_tissue_parameters_in_xcat( void );

}// namespace test_contgen



// tissue label mapper
namespace test_tlm
{

bool test_get_filepath_tissue_parameter_xml( void );
bool test_get_labels_array(void);
bool test_get_segmentation_dimensions( void );


bool test_assign_tissue_parameters_label_found( void );
bool test_assign_tissue_parameters_label_not_found( void );

bool test_map_labels_to_tissue_from_xml( void );
bool test_replace_petmr_tissue_parameters( void );

}// namespace test_tlm





