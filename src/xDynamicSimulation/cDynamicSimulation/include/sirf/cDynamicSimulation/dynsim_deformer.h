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



/*!
	\brief Utility class to deform image data with motion fields during a simulation.
	This class contains an affine offset transformation and can deform contrast generators.
*/
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
	void add_offset_deformation(const std::vector<sirf::NiftiImageData3DDeformation<float> > &vec)
	{
		displacement_offset_ = vec;
	}

	/*!
	\brief Function to deform the image in reference motion state held by a contrast generator into the motion states defined by the displacement fields.
	*/
	void deform_contrast_generator(MRContrastGenerator& mr_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);
	void deform_contrast_generator(PETContrastGenerator& pet_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> >& vec_displacement_fields);
	
	/*!
	\brief Function to resample an image to a template geometry. 
	This function is used to resample the segmentation volume to the acquisition template geometry.
	During the resampling the first the displacement offset, then the affine offset_ transformation is applied, and then the volume is resampled into the template geometry.
	*/
	sirf::NiftiImageData3D<float> resample_to_template(sirf::NiftiImageData3D<float> img, bool const use_nearest_neighbor=false) const;
	
	void set_template_rawdata(const sirf::MRAcquisitionData& ad)
	{
		sptr_mr_template_img_ = std::shared_ptr<sirf::GadgetronImagesVector>(new sirf::GadgetronImagesVector(ad));
		mr_template_available_ = true;
	}

	

protected:


	sirf::AffineTransformation<float> offset_;
	std::vector<sirf::NiftiImageData3DDeformation<float> > displacement_offset_;

	std::shared_ptr<sirf::GadgetronImagesVector> sptr_mr_template_img_;
	bool mr_template_available_ = false;

	static const std::string temp_folder_name_;

	static void deform_pet_image( sirf::STIRImageData& img, std::vector<sirf::NiftiImageData3DDeformation<float> > &vec_displacement_fields);
	
};
