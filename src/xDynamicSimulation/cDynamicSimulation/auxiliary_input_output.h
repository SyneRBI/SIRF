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
#include <gadgetron/ImageIOAnalyze.h>


#include "gadgetron_data_containers.h"

#define SHARED_FOLDER_PATH "/media/sf_SharedFolder/CCPPETMR/"


#define ANALYZE_OUTPUT_TESTPATH SHARED_FOLDER_PATH "analyze_test_output"



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


	template <typename T>
	void write_ISMRMRD_Image_to_Analyze(std::string const output_name_without_ext, ISMRMRD::Image<T> img)
	{
		std::cout << "Started writing " << output_name_without_ext << std::endl;	

		std::vector < size_t > img_dims;
		img_dims.push_back( img.getMatrixSizeX() );
		img_dims.push_back( img.getMatrixSizeY() );
		img_dims.push_back( img.getMatrixSizeZ() );
		img_dims.push_back( img.getNumberOfChannels() );

		float const pix_size_X = img.getFieldOfViewX() / (float)img_dims[0];
		float const pix_size_Y = img.getFieldOfViewY() / (float)img_dims[1];
		float const pix_size_Z = img.getFieldOfViewZ() / (float)img_dims[2];
		float const pix_size_Vec = 1.f;


		Gadgetron::hoNDArray< T > data_to_be_written( img_dims );
		size_t const num_elements = img.getNumberOfDataElements();
		
		for( size_t i=0; i<num_elements; i++)
			*(data_to_be_written.begin() + i) = *(img.begin() + i);
	

		Gadgetron::ImageIOAnalyze analyze_io( pix_size_X, pix_size_Y, pix_size_Z, pix_size_Vec);

		analyze_io.export_array(data_to_be_written, output_name_without_ext);

		std::cout << "Finished writing " << std::endl;	
	};

}

namespace mr_io{

	ISMRMRD::IsmrmrdHeader read_ismrmrd_header( std::string path_ismrmrd_h5_file_with_ext);
	sirf::AcquisitionsVector read_ismrmrd_acquisitions( std::string path_ismrmrd_h5_file_with_ext);

}