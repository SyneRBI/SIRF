/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <string>
#include <vector>

#include <ismrmrd/ismrmrd.h>

#include "sirf/cDynamicSimulation/contrastgenerator.h"

#include "sirf/cDynamicSimulation/auxiliary_input_output.h" // this header (rather the Gadgetron Base IO including Nifti) must not be included after the SIRFImageData.h headers. DONT put it into the cpp!

#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/STIR/stir_data_containers.h"

class DynamicSimulationDeformer
{

public:

	
	static void deform_contrast_generator(MRContrastGenerator& mr_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);

	static void deform_contrast_generator(PETContrastGenerator& pet_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);


	static ISMRMRD::Image< float > extract_real_part( ISMRMRD::Image< complex_float_t >& img );
	static ISMRMRD::Image< float > extract_imaginary_part( ISMRMRD::Image< complex_float_t >& img );

protected:

	static const std::string temp_folder_name_;

	static void deform_ismrmrd_image(ISMRMRD::Image< float >& img, std::vector<sirf::NiftiImageData3DDeformation<float> > &vec_displacement_fields);
	static void deform_pet_image( sirf::STIRImageData& img, std::vector<sirf::NiftiImageData3DDeformation<float> > &vec_displacement_fields);
	
	static ISMRMRD::Image< float > extract_complex_subpart( ISMRMRD::Image< complex_float_t >& img, bool const extract_real_part );

};
