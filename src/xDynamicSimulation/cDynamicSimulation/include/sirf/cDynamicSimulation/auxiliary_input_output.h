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
// #include <gadgetron/ImageIOAnalyze.h>

#include "sirf/STIR/stir_data_containers.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
 

#include "sirf/common/GeometricalInfo.h"
#include "sirf/common/ImageData.h"

#include "sirf/Reg/NiftiImageData3D.h"
	

#include "test_input_filenames.h"

#define PIX_SIZE_X 3.20f
#define PIX_SIZE_Y 3.20f
#define PIX_SIZE_Z 1.92f

// #define PIX_SIZE_X 2.0f
// #define PIX_SIZE_Y 2.0f
// #define PIX_SIZE_Z 2.0f


typedef float TimeAxisType;
typedef float SignalAxisType;
typedef std::pair<TimeAxisType, SignalAxisType> SignalPoint;
typedef std::vector< SignalPoint > SignalContainer;


namespace print_io{

	void print_voxelized_geometrical_info( const sirf::ImageData& im );
}


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
	void write_ISMRMRD_Image_to_nii(std::string const output_name_without_ext, ISMRMRD::Image<T> &img)
	{
		std::cout << "Started writing " << output_name_without_ext << std::endl;	
		std::cout << "Only first channel will be written to file. " << output_name_without_ext << std::endl;	

		// determine data size. only take first channel
		sirf::VoxelisedGeometricalInfo3D::Size data_size{ (unsigned int)img.getMatrixSizeX() , (unsigned int)img.getMatrixSizeY() , (unsigned int)img.getMatrixSizeZ() };
		
		// set spacing
		float const pix_size_X = (float)img.getFieldOfViewX() / (float)data_size[0];
		float const pix_size_Y = (float)img.getFieldOfViewY() / (float)data_size[1];
		float const pix_size_Z = (float)img.getFieldOfViewZ() / (float)data_size[2];

		sirf::VoxelisedGeometricalInfo3D::Spacing spacing{pix_size_X, pix_size_Y,  pix_size_Z};

		
		sirf::VoxelisedGeometricalInfo3D::Offset offset;
		for( int i=0; i<3; i++)
			offset[i] = (-0.5f * data_size[i] + 0.5f) * spacing[i];


		sirf::VoxelisedGeometricalInfo3D::DirectionMatrix dir_mat;
	
		dir_mat[0] = {1,0,0};
		dir_mat[1] = {0,1,0};
		dir_mat[2] = {0,0,1};

		sirf::VoxelisedGeometricalInfo3D geometry_info( offset, spacing, data_size, dir_mat);
			
		
		size_t const num_relevant_elements = data_size[0]*data_size[1]*data_size[2] ;
		
		
		std::vector< float > data( num_relevant_elements );


		for( size_t i=0; i<num_relevant_elements; i++)
		 	data[i] = std::abs(*(img.begin() + i));
	
		sirf::NiftiImageData3D<float> nifti_img(&data[0], geometry_info);
		nifti_img.write( output_name_without_ext + ".nii");

		std::cout << "Finished writing "  << output_name_without_ext << std::endl;		
	};

	template <typename T>
	void write_MVF_from_ISMRMRD_Image_to_Analyze(std::string const output_name_without_ext, ISMRMRD::Image<T> mvf)
	{

		if(mvf.getNumberOfChannels() != 3)
			throw std::runtime_error("Please pass a 3d vector field to write out.");

		size_t const Nx = mvf.getMatrixSizeX();
		size_t const Ny = mvf.getMatrixSizeY();
		size_t const Nz = mvf.getMatrixSizeZ();
		
		std::cout << "Rescaling motion fields to milimeters." <<std::endl;

		for(size_t nz=0; nz<Nz; nz++)
		for(size_t ny=0; ny<Ny; ny++)
		for(size_t nx=0; nx<Nx; nx++)
		{
			mvf(nx, ny, nz, 0) *= PIX_SIZE_Z;
			mvf(nx, ny, nz, 1) *= PIX_SIZE_Y;
			mvf(nx, ny, nz, 2) *= PIX_SIZE_X;
		}

		write_ISMRMRD_Image_to_nii<T> (output_name_without_ext, mvf);
	
	}

	template <typename T>
	void write_MVF_from_ISMRMRD_Image_to_Analyze_In_PET_Geometry(std::string const output_name_without_ext, ISMRMRD::Image<T> mvf)
	{

		if(mvf.getNumberOfChannels() != 3)
			throw std::runtime_error("Please pass a 3d vector field to write out.");

		size_t const Nx = mvf.getMatrixSizeX();
		size_t const Ny = mvf.getMatrixSizeY();
		size_t const Nz = mvf.getMatrixSizeZ();
		
		std::cout << "Rescaling motion fields to milimeters." <<std::endl;

		for(size_t nz=0; nz<Nz; nz++)
		for(size_t ny=0; ny<Ny; ny++)
		for(size_t nx=0; nx<Nx; nx++)
		{
			mvf(nx, ny, nz, 0) *= -PIX_SIZE_X;
			mvf(nx, ny, nz, 1) *= -PIX_SIZE_Y;
			mvf(nx, ny, nz, 2) *= PIX_SIZE_Z;
		}

		write_ISMRMRD_Image_to_nii<T> (output_name_without_ext, mvf);
	
	}




	void write_PET_image_to_hv( const std::string& filename_without_ext,const sirf::STIRImageData& img);


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
		      while(myfile >> val)
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