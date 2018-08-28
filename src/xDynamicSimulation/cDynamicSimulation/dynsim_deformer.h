/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <string>

#include <ismrmrd/ismrmrd.h>

#include "contrastgenerator.h"

#include "auxiliary_input_output.h" // this header (rather the Gadgetron Base IO including Nifti) must not be included after the SIRFImageData.h headers. DONT put it into the cpp!

#include "SIRFImageDataDeformation.h"


class DynamicSimulationDeformer
{

public:

	static void deform_contrast_generator(MRContrastGenerator& mr_cont_gen, SIRFImageDataDeformation& displacement_field);
	
	static ISMRMRD::Image< float > extract_real_part( ISMRMRD::Image< complex_float_t >& img );
	static ISMRMRD::Image< float > extract_imaginary_part( ISMRMRD::Image< complex_float_t >& img );

protected:

	static const std::string temp_folder_name_;

	static ISMRMRD::Image< float > extract_complex_subpart( ISMRMRD::Image< complex_float_t >& img, bool const extract_real_part );

};