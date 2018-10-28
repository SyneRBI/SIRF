/* ================================================

Author: Johannes Mayer
Date: 2018.04.05
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once


#ifdef _SIRFIMAGEDATA_H_
	std::cout << "###############################################################" <<std::endl;
	std::cout << "######################### WARNING #############################" <<std::endl;
	std::cout << "THE HEADER FILE SIRFImageData.h WAS INCLUDED ALREADY." <<std::endl;
	std::cout << "COULD LEAD TO PROBLEMS." <<std::endl;
	std::cout << "###############################################################" <<std::endl;
#endif


#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include <ismrmrd/xml.h>
#include <gadgetron/ImageIOAnalyze.h>


#include "gadgetron_data_containers.h"

#define SHARED_FOLDER_PATH "/media/sf_SharedFolder/CCPPETMR/"


#define ANALYZE_OUTPUT_TESTPATH SHARED_FOLDER_PATH "analyze_test_output"

typedef float TimeAxisType;
typedef float SignalAxisType;
typedef std::pair<TimeAxisType, SignalAxisType> SignalPoint;
typedef std::vector< SignalPoint > SignalContainer;

namespace data_io{

	template <typename T>
	void write_raw(std::string const output_name_without_ext, T* ptr_to_data, size_t const num_elements)
	{	
		std::cout<< "Writing file " << output_name_without_ext << std::endl;
		std::stringstream name_stream;
		name_stream << output_name_without_ext << ".bin";

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
	void write_ISMRMRD_Image_to_Analyze(std::string const output_name_without_ext, ISMRMRD::Image<T> &img)
	{
		std::cout << "Started writing " << output_name_without_ext << std::endl;	

		std::vector < size_t > img_dims;
		img_dims.push_back( img.getMatrixSizeX() );
		img_dims.push_back( img.getMatrixSizeY() );
		img_dims.push_back( img.getMatrixSizeZ() );
		img_dims.push_back( 1 );
		img_dims.push_back( img.getNumberOfChannels() );

		float const pix_size_X = 2.f;//img.getFieldOfViewX() / (float)img_dims[0];
		float const pix_size_Y = 2.f;//img.getFieldOfViewY() / (float)img_dims[1];
		float const pix_size_Z = 2.f;//img.getFieldOfViewZ() / (float)img_dims[2];
		// float const pix_size_X = (float)img.getFieldOfViewX() / (float)img_dims[0];
		// float const pix_size_Y = (float)img.getFieldOfViewY() / (float)img_dims[1];
		// float const pix_size_Z = (float)img.getFieldOfViewZ() / (float)img_dims[2];
		// std::cout << "########################" << std::endl;
		// std::cout << img.getFieldOfViewX() << std::endl;
		// std::cout << img.getFieldOfViewY() << std::endl;
		// std::cout << img.getFieldOfViewZ() << std::endl;
		// std::cout << "########################" << std::endl;
		float const pix_size_U = 0.f;
		float const pix_size_Vec = 1.f;


		Gadgetron::hoNDArray< T > data_to_be_written( img_dims );
		size_t const num_elements = img.getNumberOfDataElements();
		
		for( size_t i=0; i<num_elements; i++)
			*(data_to_be_written.begin() + i) = *(img.begin() + i);
	

		Gadgetron::ImageIOAnalyze analyze_io( pix_size_X, pix_size_Y, pix_size_Z, pix_size_U, pix_size_Vec);

		analyze_io.export_array(data_to_be_written, output_name_without_ext);

		std::cout << "Finished writing "  << output_name_without_ext << std::endl;		
	};

	template <typename T>
	std::vector< T > read_single_column_txt( const std::string& filename_txt_without_ext )
	{
		std::string const filename_with_ext = filename_txt_without_ext + ".txt";
		std::vector <T> output;
		std::string line; 
		std::ifstream myfile (filename_with_ext);
		if (myfile.is_open())
		  {
		  	while ( getline (myfile,line) )
		    {
		      T val;
		      myfile >> val;
		      output.push_back(val);
		    }
		    myfile.close();
		  }

		 else
			 throw std::runtime_error( "Unable to open file "+filename_with_ext);

		return output;
	}

	SignalContainer read_surrogate_signal( const std::string& filename_time_axis, const std::string& filename_signal );

}

namespace mr_io{

	ISMRMRD::IsmrmrdHeader read_ismrmrd_header( std::string path_ismrmrd_h5_file_with_ext);
	sirf::AcquisitionsVector read_ismrmrd_acquisitions( std::string path_ismrmrd_h5_file_with_ext);

}