/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 - 2022 Physikalisch-Technische Bundesanstalt (PTB)

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

/*!
\file
\ingroup Simulation
\brief Classes and utilities for deforming MR and PET images during simulation.

\author Johannes Mayer
\author SyneRBI
*/

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


namespace sirf{



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

} //namespace sirf