/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>
#include <stdio.h>
#include <iostream>

#include <ismrmrd/ismrmrd.h>

#include "phantom_input.h"

bool test_read_h5_segmentation_correct_dims( std::string h5_filename_with_suffix);
bool test_read_h5_segmentation_correct_content( std::string h5_filename_with_suffix);
