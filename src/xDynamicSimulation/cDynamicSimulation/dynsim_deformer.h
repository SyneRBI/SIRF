/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <string>

#include "contrastgenerator.h"

#include "auxiliary_input_output.h" // this header (rather the Gadgetron Base IO including Nifti) must not be included after the SIRFImageData.h headers. DONT put it into the cpp!

#include "SIRFImageDataDeformation.h"


class DynamicSimulationDeformer
{

public:

	static void deform_contrast_generator(MRContrastGenerator& mr_cont_gen, SIRFImageDataDeformation& displacement_field);

protected:

	static const std::string temp_folder_name_;

};