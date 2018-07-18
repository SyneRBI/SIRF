/* ================================================

Author: Johannes Mayer
Date: 2018.04.05
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once


#include <string>
#include <sstream>
#include <fstream>

#include <ismrmrd/xml.h>

#include "gadgetron_data_containers.h"

#define SHARED_FOLDER_PATH std::string("/media/sf_SharedFolder/CCPPETMR/")

namespace data_io{

template <typename T>
	void write_raw(std::string const output_name_without_ext, T* ptr_to_data, size_t const num_elements)
	{	
		std::cout<< "Writing file " << output_name_without_ext << std::endl;
		std::stringstream name_stream;
		name_stream << output_name_without_ext << ".raw";

		std::vector <T> buffer;
		buffer.resize(num_elements);

		for( size_t i=0; i<num_elements; i++)
			buffer[i] = *(ptr_to_data + i);
		

		std::ofstream out;
		out.open( name_stream.str().c_str(), std::ios::out | std::ios::binary);

		out.write( reinterpret_cast<char*> (buffer.data()), buffer.size()*sizeof(T));

		out.close();

		std::cout<< "Finished writing file " << name_stream.str() << std::endl;

	};

}

namespace mr_io{

	ISMRMRD::IsmrmrdHeader read_ismrmrd_header( std::string path_ismrmrd_h5_file_with_ext);
	AcquisitionsVector read_ismrmrd_acquisitions( std::string path_ismrmrd_h5_file_with_ext);
}