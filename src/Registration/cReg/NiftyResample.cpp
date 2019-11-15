/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2019 University College London

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

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
\ingroup Registration
\brief Resampling class based on nifty resample

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftyResample.h"
#include "sirf/Reg/NiftiImageData3DTensor.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/NiftyMoMo/BSplineTransformation.h"
#include <_reg_resampling.h>
#include <_reg_globalTrans.h>
#include <_reg_tools.h>
#include <memory>

using namespace sirf;

template<class dataType>
void NiftyResample<dataType>::process()
{
    std::cout << "\n\nStarting resampling...\n\n";

    // Check that all the required information has been entered
    this->check_parameters();

    // Get reference and floating images as NiftiImageData
    set_up_input_images();

    // Setup output image
    set_up_output_image();
    // If no transformations, use identity.
    if (this->_transformations.size() == 0) {
        std::cout << "\nNo transformations set, using identity.\n";
        this->_transformations.push_back(std::make_shared<AffineTransformation<float> >());
    }

    // If there are multiple transformations, compose them into single transformation.
    NiftiImageData3DDeformation<dataType> transformation =
            NiftiImageData3DDeformation<dataType>::compose_single_deformation(this->_transformations, *this->_reference_image_nifti_sptr);

    // If we're not doing the adjoint
    if (!this->_do_adjoint)
        this->transformation(transformation);

    // If we are doing the adjoint
    else
        this->transformation_adjoint(transformation);

    // The output should be a clone of the reference image, with data filled in from the nifti image
    this->_output_image_sptr = this->_reference_image_sptr->clone();
    this->_output_image_sptr->fill(*this->_output_image_nifti_sptr);

    std::cout << "\n\nResampling finished!\n\n";
}

template<class dataType>
void NiftyResample<dataType>::set_up_input_images()
{
    // Try to dynamic cast from ImageData to NiftiImageData. This will only succeed if original type was NiftiImageData
    this->_reference_image_nifti_sptr = std::dynamic_pointer_cast<const NiftiImageData<dataType> >(this->_reference_image_sptr);
    this->_floating_image_nifti_sptr  = std::dynamic_pointer_cast<const NiftiImageData<dataType> >(this->_floating_image_sptr);

    // If either is a null pointer, it means that a different image type was supplied (e.g., STIRImageData).
    // In this case, construct a NiftiImageData
    if (!this->_reference_image_nifti_sptr)
        this->_reference_image_nifti_sptr = std::make_shared<const NiftiImageData<dataType> >(*this->_reference_image_sptr);
    if (!this->_floating_image_nifti_sptr)
        this->_floating_image_nifti_sptr = std::make_shared<const NiftiImageData<dataType> >(*this->_floating_image_sptr);
}

template<class dataType>
void NiftyResample<dataType>::set_up_output_image()
{
    // The output is a mixture between the reference and floating images.
    const nifti_image * const ref_ptr = this->_reference_image_nifti_sptr->get_raw_nifti_sptr().get();
    const nifti_image * const flo_ptr = this->_floating_image_nifti_sptr->get_raw_nifti_sptr().get();

    // Start creating new output as the header from the reference image
    nifti_image * output_ptr = nifti_copy_nim_info(ref_ptr);

    // Put in the required info from the floating image
    output_ptr->cal_min     = flo_ptr->cal_min;
    output_ptr->cal_max     = flo_ptr->cal_max;
    output_ptr->scl_slope   = flo_ptr->scl_slope;
    output_ptr->scl_inter   = flo_ptr->scl_inter;
    output_ptr->datatype    = flo_ptr->datatype;
    output_ptr->intent_code = flo_ptr->intent_code;
    output_ptr->intent_p1   = flo_ptr->intent_p1;
    output_ptr->intent_p2   = flo_ptr->intent_p2;
    output_ptr->datatype    = flo_ptr->datatype;
    output_ptr->nbyper      = flo_ptr->nbyper;
    memset(output_ptr->intent_name, 0, 16);
    strcpy(output_ptr->intent_name,flo_ptr->intent_name);
    output_ptr->nvox = this->_reference_image_nifti_sptr->get_num_voxels();

    // Allocate the data
    output_ptr->data = static_cast<void *>(calloc(output_ptr->nvox, unsigned(output_ptr->nbyper)));

    // Create NiftiImageData from nifti_image
    this->_output_image_nifti_sptr = std::make_shared<NiftiImageData<dataType> >(*output_ptr);
}

template<class dataType>
void NiftyResample<dataType>::transformation(NiftiImageData3DDeformation<dataType> &transformation)
{
    // Annoyingly NiftyReg doesn't mark floating image as const, so need to copy (could do a naughty C-style cast?)
    NiftiImageData<dataType> flo = *this->_floating_image_nifti_sptr;

    reg_resampleImage(flo.get_raw_nifti_sptr().get(),
                      this->_output_image_nifti_sptr->get_raw_nifti_sptr().get(),
                      transformation.get_raw_nifti_sptr().get(),
                      NULL,
                      this->_interpolation_type,
                      this->_padding_value);
}

template<class dataType>
void NiftyResample<dataType>::transformation_adjoint(NiftiImageData3DDeformation<dataType> &transformation)
{
    // SINC currently not supported in NiftyMoMo
    if (this->_interpolation_type == Resample<dataType>::SINC)
        throw std::runtime_error("NiftyMoMo does not currently support SINC interpolation");

    // Get deformation spacing
    float control_point_grid_spacing[3] = {
        transformation.get_raw_nifti_sptr()->dx,
        transformation.get_raw_nifti_sptr()->dy,
        transformation.get_raw_nifti_sptr()->dz
    };

    // Annoyingly NiftyReg doesn't mark reference image as const, so need to copy (could do a naughty C-style cast?)
    NiftiImageData<dataType> ref = *this->_reference_image_nifti_sptr;
    nifti_image *ref_ptr = ref.get_raw_nifti_sptr().get();

    // Output image is a copy of the floating image
    this->_output_image_nifti_sptr = this->_floating_image_nifti_sptr->clone();

    // Get the raw data of the transformation
    nifti_image * def_ptr = transformation.get_raw_nifti_sptr().get();

    // Need some weights. Same as original image but filled with 1's
    NiftiImageData<dataType> ref_weights = ref;
    ref_weights.fill(1.f);
    NiftiImageData<dataType> warped_weights = *this->_output_image_nifti_sptr;
    warped_weights.fill(1.f);

    NiftyMoMo::BSplineTransformation adjoint_transformation(
                ref.get_raw_nifti_sptr().get(),
                /*num levels to perform*/1U,
                control_point_grid_spacing);

    adjoint_transformation.set_interpolation(this->_interpolation_type);
    adjoint_transformation.SetParameters(static_cast<dataType*>(def_ptr->data), false);
    adjoint_transformation.TransformImageAdjoint(ref_ptr,
                                                 ref_weights.get_raw_nifti_sptr().get(),
                                                 this->_output_image_nifti_sptr->get_raw_nifti_sptr().get(),
                                                 warped_weights.get_raw_nifti_sptr().get());
}

namespace sirf {
template class NiftyResample<float>;
}

