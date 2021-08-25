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
#include "sirf/Gadgetron/gadgetron_data_containers.h"

class DynamicSimulationDeformer
{

public:

	void deform_contrast_generator(MRContrastGenerator& mr_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);
	void deform_contrast_generator(PETContrastGenerator& pet_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);

	void set_template_rawdata(const sirf::MRAcquisitionData& ad)
	{
		sptr_mr_template_img_ = std::shared_ptr<sirf::GadgetronImageData>(new sirf::GadgetronImagesVector(ad));
		mr_template_available_ = true;
	}

protected:

	std::shared_ptr<sirf::GadgetronImageData> sptr_mr_template_img_;
	bool mr_template_available_ = false;

	static const std::string temp_folder_name_;

	static void deform_pet_image( sirf::STIRImageData& img, std::vector<sirf::NiftiImageData3DDeformation<float> > &vec_displacement_fields);

};
