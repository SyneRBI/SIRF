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

#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"

#include "sirf/STIR/stir_data_containers.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"

class DynamicSimulationDeformer
{

public:

	DynamicSimulationDeformer(){

		std::array<float, 3> id_translation{0,0,0};
		std::array<float, 3> id_euler{0,0,0};

		offset_ = sirf::AffineTransformation<float>(id_translation, id_euler);

	}

	void set_offset_transformation(const sirf::AffineTransformation<float>& trafo)
	{
		offset_ = trafo;
	}

	void deform_contrast_generator(MRContrastGenerator& mr_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);
	void deform_contrast_generator(PETContrastGenerator& pet_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);
	
	sirf::NiftiImageData3D<float> resample_to_template(const sirf::NiftiImageData3D<float>& img) const;
	
	void set_template_rawdata(const sirf::MRAcquisitionData& ad)
	{
		sptr_mr_template_img_ = std::shared_ptr<sirf::GadgetronImagesVector>(new sirf::GadgetronImagesVector(ad));
		mr_template_available_ = true;
	}

protected:

	sirf::AffineTransformation<float> offset_;

	std::shared_ptr<sirf::GadgetronImagesVector> sptr_mr_template_img_;
	bool mr_template_available_ = false;

	static const std::string temp_folder_name_;

	static void deform_pet_image( sirf::STIRImageData& img, std::vector<sirf::NiftiImageData3DDeformation<float> > &vec_displacement_fields);
	
	std::shared_ptr<sirf::AffineTransformation<float> > compute_shift_to_center(const sirf::NiftiImageData3D<float>& img) 
	{
	
		int const num_dims = 7;
		std::array<float, num_dims> spacings;
		std::array<int, num_dims> dimensions;
		
		for (int i=0; i<= num_dims; ++i)
		{
			dimensions[i] = img.get_raw_nifti_sptr()->dim[i];
            spacings[i] = img.get_raw_nifti_sptr()->pixdim[i];
		}

		std::array<float, 3> id_euler{0,0,0};
		std::array<float, 3> shift_translation{0,0,0};
		
		shift_translation[2] = -1 * dimensions[3]/2.f * spacings[3];
		sirf::AffineTransformation<float> offset(shift_translation, id_euler);
		offset_ = offset;

		return std::shared_ptr<sirf::AffineTransformation<float> >(new sirf::AffineTransformation<float>(offset));
 	}
};
