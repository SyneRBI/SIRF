/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#define H5_TEST_PATH "Testing/TestData/h5_testfile_cube_size3.h5"



#include "tests_tissueparameters.h"
#include "tests_contrastgenerator.h"
#include "tests_phantom_input.h"

// Function declarations in order to collect tests for module dynamicssimulation

void run_tests_tissueparameters( void );

void run_tests_phantom_input( void );

void run_tests_contrastgenerator( void );
