/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "dynsim_deformer.h"

#include <stdexcept>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>

#include "SIRFRegNiftyResample.h"

std::string const DynamicSimulationDeformer::temp_folder_name_ = "/tmp/tmp_img_data_to_deform";

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
			data_io::write_ISMRMRD_Image_to_Analyze< complex_float_t > (filename_temp_img, vect_img_data[i_cont]);

			filename_temp_img += ".hdr";

		    SIRFImageData img_to_deform(filename_temp_img);

		    SIRFRegNiftyResample resampler; 

		    resampler.set_interpolation_type_to_linear();
    		resampler.set_reference_image(img_to_deform);
    		resampler.set_floating_image(img_to_deform);

			resampler.set_displacement_field(displacement_field);

			resampler.update();

		}
	}
	else
	{
		throw std::runtime_error("The temporary folder creation to store image data to deform failed. Please choose a path which does not exist yet or to which you have access.");
	}

	bool  dummy_variable = boost::filesystem::remove_all(temp_folder_name_);
}

