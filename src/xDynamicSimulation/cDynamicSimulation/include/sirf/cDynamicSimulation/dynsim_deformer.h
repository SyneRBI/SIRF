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
#include "sirf/Reg/NiftiImageData3DDisplacement.h"

#include "sirf/STIR/stir_data_containers.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"

class DynamicSimulationDeformer
{

public:

	void deform_contrast_generator(MRContrastGenerator& mr_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);
	void deform_contrast_generator(PETContrastGenerator& pet_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);

	void set_template_rawdata(const sirf::MRAcquisitionData& ad)
	{
		sptr_mr_template_img_ = std::shared_ptr<sirf::GadgetronImagesVector>(new sirf::GadgetronImagesVector(ad));
		mr_template_available_ = true;
	}

protected:

	std::shared_ptr<sirf::GadgetronImagesVector> sptr_mr_template_img_;
	bool mr_template_available_ = false;

	static const std::string temp_folder_name_;

	static void deform_pet_image( sirf::STIRImageData& img, std::vector<sirf::NiftiImageData3DDeformation<float> > &vec_displacement_fields);

	std::shared_ptr<sirf::NiftiImageData3DDeformation<float> > compute_shift_to_center(const sirf::NiftiImageData3D<float>& img) const
	{
		sirf::NiftiImageData3D<float> shift_x(img), shift_y(img), shift_z(img);		
		
		int const num_dims = 7;
		std::array<float, num_dims> spacings;
		std::array<int, num_dims> dimensions;
		
		for (int i=0; i<= num_dims; ++i)
		{
			dimensions[i] = shift_z.get_raw_nifti_sptr()->dim[i];
            spacings[i] = shift_z.get_raw_nifti_sptr()->pixdim[i];
		}

		shift_x.fill(0.f);	
		shift_y.fill(0.f);	
		shift_z.fill(-1 * dimensions[3] * spacings[3]);
		
		sirf::NiftiImageData3DDisplacement<float> shift_vf(shift_x,shift_y,shift_z);
		return std::make_shared<sirf::NiftiImageData3DDeformation<float> >(shift_vf);

	}
};
