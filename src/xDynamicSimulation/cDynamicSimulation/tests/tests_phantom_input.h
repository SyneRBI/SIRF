 /* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>
#include <stdio.h>
#include <iostream>


bool test_read_1D_dataset_from_h5( std::string h5_filename_with_suffix );
bool test_read_geometrical_info_from_h5( std::string h5_filename_with_suffix );
bool test_read_segmentation_to_nifti( std::string h5_filename_with_suffix );
bool test_read_motionfield_to_nifti(  std::string h5_filename_with_suffix );
