/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <ismrmrd/ismrmrd.h>

ISMRMRD::NDArray< unsigned int > read_segmentation_from_h5( std::string const h5_filename_with_suffix);