/* ================================================

Author: Johannes Mayer
Date: 2018.04.05
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "auxiliary_input_output.h"


#include <ismrmrd/ismrmrd.xml.h>



// ++++++++++++++++++ mr_IO ++++++++++++++++++

ISMRMRD::IsmrmrdHeader mr_io::read_ismrmrd_header( std::string path_ismrmrd_h5_file_with_ext)
{

	ISMRMRD::Dataset d(path_ismrmrd_h5_file_with_ext.c_str(),"dataset", false);

	std::string xml;
	d.readHeader(xml);

	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(), hdr);
	
	return hdr;

}


// ++++++++++++++++++ pet_IO ++++++++++++++++++