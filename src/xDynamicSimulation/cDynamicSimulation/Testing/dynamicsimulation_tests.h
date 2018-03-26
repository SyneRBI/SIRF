/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once



#define XML_TEST_PATH "Testing/TestData/test_TissueParameters_XML.xml" 
#define H5_TEST_PATH "Testing/TestData/h5_testfile_cube_size3.h5"



#include "tests_tissueparameters.h"
#include "tests_h5_reader.h"

// Function declarations in order to collect tests for module dynamicssimulation

void run_tests_tissueparameters( void );

void run_tests_h5_reader( void );
