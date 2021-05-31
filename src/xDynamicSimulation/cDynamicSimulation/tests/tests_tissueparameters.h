/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

This file collects all testjobs relevant for the tissue
parameter structs and the xml parser filling them.

================================================ */

#pragma once


#include <string>
#include <stdio.h>
#include <iostream>


#include "sirf/cDynamicSimulation/tissueparameters.h"



bool test_allocate_MRTissueParameter_successful(void);
bool test_allocate_PETTissueParameter_successful(void);
bool test_allocate_TissueParameter_successful(void);

bool test_TissueParameter_algebra( void );

bool test_get_MRTissueParameter_from_ptree(void);
bool test_get_PETTissueParameter_from_ptree(void);

bool test_exception_throw_if_node_not_exists(void);

bool test_read_TissueParameter_label_from_xml( std::string const xml_filepath );

//TissueParameterList get_test_tissue_parameter_list( void );


bool test_check_label_uniqueness_fails();
bool test_check_label_uniqueness_true();

