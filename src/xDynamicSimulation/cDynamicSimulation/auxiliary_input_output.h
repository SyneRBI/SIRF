/* ================================================

Author: Johannes Mayer
Date: 2018.04.05
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once


#include <string>


#include <ismrmrd/xml.h>



namespace mr_io{

	ISMRMRD::IsmrmrdHeader read_ismrmrd_header( std::string path_ismrmrd_h5_file_with_ext);

}