/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "dynsim_deformer.h"


#include <stdexcept>
#include <sstream>
#include <complex>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>

#include "SIRFRegNiftyResample.h"

std::string const DynamicSimulationDeformer::temp_folder_name_ = "/tmp/tmp_img_data_to_deform";



ISMRMRD::Image< float > DynamicSimulationDeformer::extract_real_part(  ISMRMRD::Image< complex_float_t >& complex_img )
{
	return DynamicSimulationDeformer::extract_complex_subpart( complex_img, true);
}


ISMRMRD::Image< float > DynamicSimulationDeformer::extract_imaginary_part(  ISMRMRD::Image< complex_float_t >& complex_img )
{
	return DynamicSimulationDeformer::extract_complex_subpart( complex_img, false);
}

ISMRMRD::Image< float > DynamicSimulationDeformer::extract_complex_subpart( ISMRMRD::Image< complex_float_t >& complex_img, bool const extract_real_part )
{

	auto cplx_img_header = complex_img.getHead();
	cplx_img_header.data_type = 5; 

	ISMRMRD::Image< float > sub_img;
	sub_img.setHead(cplx_img_header);
	sub_img.resize( sub_img.getMatrixSizeX(), sub_img.getMatrixSizeY(), sub_img.getMatrixSizeZ(), sub_img.getNumberOfChannels() );


	for( size_t i=0; i<complex_img.getNumberOfDataElements(); i++)
	{
		if (extract_real_part)
			*(sub_img.begin() + i) = std::real(  *(complex_img.begin() + i ) );
		else
			*(sub_img.begin() + i) = std::imag(  *(complex_img.begin() + i ) );
	}

	return sub_img;
}



void DynamicSimulationDeformer::deform_contrast_generator(MRContrastGenerator& mr_cont_gen, SIRFImageDataDeformation& displacement_field)
{

	boost::filesystem::path temp_dir_name(temp_folder_name_);
	bool const temp_folder_creation_successful = boost::filesystem::create_directories(temp_dir_name);

	std::stringstream namestream_temp_img_output;
	namestream_temp_img_output << temp_folder_name_ << "/temp_img_data";

	if( temp_folder_creation_successful )
	{
		std::vector< ISMRMRD::Image< complex_float_t> >  vect_img_data = mr_cont_gen.get_contrast_filled_volumes();

		for(size_t i_cont=0; i_cont<vect_img_data.size(); i_cont++)
		{
			std::string filename_temp_img = namestream_temp_img_output.str();
			auto real_image_part =  DynamicSimulationDeformer::extract_real_part(vect_img_data[i_cont]);

			data_io::write_ISMRMRD_Image_to_Analyze< float > (filename_temp_img, real_image_part);

			filename_temp_img += ".hdr";

		    SIRFImageData img_to_deform(filename_temp_img);

		    SIRFRegNiftyResample resampler; 

		    resampler.set_interpolation_type_to_linear();
    		resampler.set_reference_image(img_to_deform);
    		resampler.set_floating_image(img_to_deform);

			resampler.set_displacement_field(displacement_field);

			resampler.update();

			auto deformed_img = resampler.get_output();
			auto deformed_img_as_nifti = *(deformed_img.get_image_as_nifti());
			
			if( deformed_img_as_nifti.nvox != real_image_part.getNumberOfDataElements() )
				throw std::runtime_error("Something went wrong during the resampling. The output image and input image have different number of voxels.");

			for( size_t i_vox=0; i_vox< deformed_img_as_nifti.nvox; i_vox++)			
				*(real_image_part.begin() + i_vox) = ((float*) deformed_img_as_nifti.data)[i_vox];


			data_io::write_ISMRMRD_Image_to_Analyze< float > (filename_temp_img, real_image_part);

		}
	}
	else
	{
		throw std::runtime_error("The temporary folder creation to store image data to deform failed. Please choose a path which does not exist yet or to which you have access.");
	}

	bool  dummy_variable = boost::filesystem::remove_all(temp_folder_name_);
}

