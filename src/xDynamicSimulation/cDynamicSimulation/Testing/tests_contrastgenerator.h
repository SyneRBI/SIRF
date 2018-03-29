/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once


#include <string>
#include <stdio.h>
#include <iostream>


#include "auxiliary_testing_functions.h"

#include "tissuelabelmapper.h"
#include "tissueparameters.h"
#include "contrastgenerator.h"


namespace test_contgen
{

bool test_mr_constructor( void );
bool test_mr_set_get_rawdata_header_path( void );
bool test_mr_read_rawdata_header( void );

}// namespace test_tlm



// tissue label mapper
namespace test_tlm
{

bool test_get_filepath_tissue_parameter_xml( void );
bool test_get_labels_array(void);

bool test_assign_tissue_parameters_label_found( void );
bool test_assign_tissue_parameters_label_not_found( void );

bool test_map_labels_to_tissue_from_xml( void );

}// namespace test_tlm





