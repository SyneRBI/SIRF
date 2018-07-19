/* ================================================

Author: Johannes Mayer
Date: 2018.04.05
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "auxiliary_input_output.h"


#include <ismrmrd/dataset.h>
#include <ismrmrd/ismrmrd.h>


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

AcquisitionsVector mr_io::read_ismrmrd_acquisitions( std::string path_ismrmrd_h5_file_with_ext )
{

	ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header( path_ismrmrd_h5_file_with_ext );
	std::ostringstream out;
	serialize(hdr, out);

	std::cout<< "Started reading acquisitions from " << path_ismrmrd_h5_file_with_ext << std::endl;

	ISMRMRD::Dataset d(path_ismrmrd_h5_file_with_ext.c_str(),"dataset", false);

	uint32_t num_acquis = d.getNumberOfAcquisitions();

	AcquisitionsVector acq_vec(out.str());
	ISMRMRD::Acquisition acq;

	for( uint32_t i_acqu=0; i_acqu<num_acquis; i_acqu++)
	{
		d.readAcquisition( i_acqu, acq);
		acq_vec.append_acquisition( acq );

	}

	std::cout<< "Finished reading acquisitions from " << path_ismrmrd_h5_file_with_ext << std::endl;

	return acq_vec;

}



// ++++++++++++++++++ pet_IO ++++++++++++++++++