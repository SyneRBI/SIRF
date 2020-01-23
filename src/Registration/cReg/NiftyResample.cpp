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
void convert_to_NiftiImageData_if_not_already(std::shared_ptr<const NiftiImageData<dataType> > &output_sptr, const std::shared_ptr<const ImageData> &input_sptr)
{
    // Try to dynamic cast from ImageData to (const) NiftiImageData. This will only succeed if original type was NiftiImageData
    output_sptr = std::dynamic_pointer_cast<const NiftiImageData<dataType> >(input_sptr);
    // If output is a null pointer, it means that a different image type was supplied (e.g., STIRImageData).
    // In this case, construct a NiftiImageData
    if (!output_sptr)
        output_sptr = std::make_shared<const NiftiImageData<dataType> >(*input_sptr);
}

template<class dataType>
void convert_to_NiftiImageData_if_not_already(std::shared_ptr<NiftiImageData<dataType> > &output_sptr, const std::shared_ptr<ImageData> &input_sptr)
{
    // Try to dynamic cast from ImageData to (const) NiftiImageData. This will only succeed if original type was NiftiImageData
    output_sptr = std::dynamic_pointer_cast<NiftiImageData<dataType> >(input_sptr);
    // If output is a null pointer, it means that a different image type was supplied (e.g., STIRImageData).
    // In this case, construct a NiftiImageData
    if (!output_sptr)
        output_sptr = std::make_shared<NiftiImageData<dataType> >(*input_sptr);
}

template<class dataType>
void NiftyResample<dataType>::set_up()
{
    // If set up has already been called, nothing to do.
    if (!this->_need_to_set_up)
        return;

    // Check that all the required information has been entered
    this->check_parameters();

    // Get reference and floating images as NiftiImageData
    set_up_input_images();

    // If no transformations, use identity.
    if (this->_transformations.size() == 0) {
        std::cout << "\nNo transformations set, using identity.\n";
        this->_transformations.push_back(std::make_shared<AffineTransformation<float> >());
    }

    // If there are multiple transformations, compose them into single transformation.
    // If forward, use the reference. If adjoint, use the floating
    this->_deformation_sptr = std::make_shared<NiftiImageData3DDeformation<dataType> >(
            NiftiImageData3DDeformation<dataType>::compose_single_deformation(
                this->_transformations,*this->_reference_image_nifti_sptr));

    this->_need_to_set_up = false;
}

template<class dataType>
void NiftyResample<dataType>::set_up_forward()
{
    // If set up has already been called, nothing to do.
    if (!this->_need_to_set_up_forward)
        return;

    // Call base level
    set_up();

    // Setup output image
    set_up_output_image(_output_image_forward_nifti_sptr, _reference_image_nifti_sptr, _floating_image_nifti_sptr);

    this->_need_to_set_up_forward = false;
}

template<class dataType>
void NiftyResample<dataType>::set_up_adjoint()
{
    // If set up has already been called, nothing to do.
    if (!this->_need_to_set_up_adjoint)
        return;

    // Call base level
    set_up();

    // SINC currently not supported in NiftyMoMo
    if (this->_interpolation_type == Resample<dataType>::SINC)
        throw std::runtime_error("NiftyMoMo does not currently support SINC interpolation");

    // Setup output image
    set_up_output_image(_output_image_adjoint_nifti_sptr, _floating_image_nifti_sptr, _reference_image_nifti_sptr);

    // Get deformation spacing
    float control_point_grid_spacing[3] = {
        _deformation_sptr->get_raw_nifti_sptr()->dx,
        _deformation_sptr->get_raw_nifti_sptr()->dy,
        _deformation_sptr->get_raw_nifti_sptr()->dz
    };

    // Get the deformation field
    nifti_image *def_ptr = _deformation_sptr->get_raw_nifti_sptr().get();
    // Copy of reference image
    NiftiImageData<dataType> ref = *this->_reference_image_nifti_sptr;
    nifti_image *ref_ptr = ref.get_raw_nifti_sptr().get();

    _adjoint_transformer_sptr = std::make_shared<NiftyMoMo::BSplineTransformation>(
                ref_ptr,
                /*num levels to perform*/1U,
                control_point_grid_spacing);

    _adjoint_transformer_sptr->set_interpolation(this->_interpolation_type);
//    _adjoint_transformer_sptr->SetParameters(static_cast<dataType*>(def_ptr->data), false);
    _adjoint_transformer_sptr->SetPaddingValue(this->_padding_value);
    _adjoint_transformer_sptr->setDVF(def_ptr);

    // Set up the input and output weights
    _adjoint_input_weights_sptr = this->_reference_image_nifti_sptr->clone();
    _adjoint_input_weights_sptr->fill(1.f);
    _adjoint_output_weights_sptr = this->_output_image_adjoint_nifti_sptr->clone();

    this->_need_to_set_up_adjoint = false;
}

template<class dataType>
void NiftyResample<dataType>::process()
{
    this->_output_image_sptr = forward(this->_floating_image_sptr);
}

template<class dataType>
void NiftyResample<dataType>::set_up_input_images()
{
    convert_to_NiftiImageData_if_not_already(this->_reference_image_nifti_sptr, this->_reference_image_sptr);
    convert_to_NiftiImageData_if_not_already(this->_floating_image_nifti_sptr,  this->_floating_image_sptr );
}

template<class dataType>
void NiftyResample<dataType>::set_up_output_image(std::shared_ptr<NiftiImageData<dataType> > &output_sptr,
                                                  const std::shared_ptr<const NiftiImageData<dataType> > im_for_shape_sptr,
                                                  const std::shared_ptr<const NiftiImageData<dataType> > im_for_metadata_sptr)
{
    // The output is a mixture between the reference and floating images.
    const nifti_image * const im_for_shape_ptr = im_for_shape_sptr->get_raw_nifti_sptr().get();
    const nifti_image * const im_for_metadata_ptr = im_for_metadata_sptr->get_raw_nifti_sptr().get();

    // Start creating new output as the header from the reference image
    nifti_image * output_ptr = nifti_copy_nim_info(im_for_shape_ptr);

    // Put in the required info from the floating image
    output_ptr->cal_min     = im_for_metadata_ptr->cal_min;
    output_ptr->cal_max     = im_for_metadata_ptr->cal_max;
    output_ptr->scl_slope   = im_for_metadata_ptr->scl_slope;
    output_ptr->scl_inter   = im_for_metadata_ptr->scl_inter;
    output_ptr->datatype    = im_for_metadata_ptr->datatype;
    output_ptr->intent_code = im_for_metadata_ptr->intent_code;
    output_ptr->intent_p1   = im_for_metadata_ptr->intent_p1;
    output_ptr->intent_p2   = im_for_metadata_ptr->intent_p2;
    output_ptr->datatype    = im_for_metadata_ptr->datatype;
    output_ptr->nbyper      = im_for_metadata_ptr->nbyper;
    memset(output_ptr->intent_name, 0, 16);
    strcpy(output_ptr->intent_name,im_for_metadata_ptr->intent_name);
    output_ptr->nvox = im_for_shape_ptr->nvox;

    // Allocate the data
    output_ptr->data = static_cast<void *>(calloc(output_ptr->nvox, unsigned(output_ptr->nbyper)));

    // Create NiftiImageData from nifti_image
    output_sptr = std::make_shared<NiftiImageData<dataType> >(*output_ptr);

    // Delete the original pointer
    nifti_image_free(output_ptr);
}

template<class dataType>
std::shared_ptr<ImageData> NiftyResample<dataType>::forward(const std::shared_ptr<const ImageData> input_sptr)
{
    // Call the set up
    set_up_forward();

    std::shared_ptr<ImageData> output_sptr = this->_reference_image_sptr->clone();
    forward(output_sptr, input_sptr);

    return output_sptr;
}

template<class dataType>
void NiftyResample<dataType>::forward(std::shared_ptr<ImageData> output_sptr, const std::shared_ptr<const ImageData> input_sptr)
{
    // Call the set up
    set_up_forward();

    // Get the input image as NiftiImageData
    // Unfortunately need to clone input image, as not marked as const in NiftyReg
    std::shared_ptr<NiftiImageData<dataType> > input_nifti_sptr, output_nifti_sptr;
    convert_to_NiftiImageData_if_not_already(input_nifti_sptr,  input_sptr->clone() );
    convert_to_NiftiImageData_if_not_already(output_nifti_sptr,  output_sptr);

    // Check that the metadata match
    if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(*input_nifti_sptr,*_floating_image_nifti_sptr,false))
        throw std::runtime_error("NiftyResample::forward: Metadata of input image should match floating image.");
    if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(*output_nifti_sptr,*_reference_image_nifti_sptr,false))
        throw std::runtime_error("NiftyResample::forward: Metadata of output image should match reference image.");

    reg_resampleImage(input_nifti_sptr->get_raw_nifti_sptr().get(),
                      this->_output_image_forward_nifti_sptr->get_raw_nifti_sptr().get(),
                      _deformation_sptr->get_raw_nifti_sptr().get(),
                      nullptr,
                      this->_interpolation_type,
                      this->_padding_value);

    // The output should be a clone of the reference image, with data filled in from the nifti image
    output_sptr->fill(*this->_output_image_forward_nifti_sptr);
    this->_output_image_sptr = output_sptr;
}

template<class dataType>
std::shared_ptr<ImageData> NiftyResample<dataType>::adjoint(const std::shared_ptr<const ImageData> input_sptr)
{
    // Call the set up
    set_up_adjoint();

    std::shared_ptr<ImageData> output_sptr = this->_floating_image_sptr->clone();
    adjoint(output_sptr, input_sptr);

    return output_sptr;
}

template<class dataType>
void NiftyResample<dataType>::adjoint(std::shared_ptr<ImageData> output_sptr, const std::shared_ptr<const ImageData> input_sptr)
{
    // Call the set up
    set_up_adjoint();

    // Get the input image as NiftiImageData
    // Unfortunately need to clone input image, as not marked as const in NiftyReg
    std::shared_ptr<NiftiImageData<dataType> > input_nifti_sptr, output_nifti_sptr;
    convert_to_NiftiImageData_if_not_already(input_nifti_sptr,  input_sptr->clone() );
    convert_to_NiftiImageData_if_not_already(output_nifti_sptr, output_sptr         );

    // Check that the metadata match
    if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(*input_nifti_sptr,*_reference_image_nifti_sptr,false))
        throw std::runtime_error("NiftyResample::adjoint: Metadata of input image should match reference image.");
    if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(*output_nifti_sptr,*_floating_image_nifti_sptr,false))
        throw std::runtime_error("NiftyResample::adjoint: Metadata of output image should match floating image.");

    // set output to zero
    this->_output_image_adjoint_nifti_sptr->fill(0.f);

    _adjoint_transformer_sptr->
            TransformImageAdjoint(input_nifti_sptr->get_raw_nifti_sptr().get(),
                                  _adjoint_input_weights_sptr->get_raw_nifti_sptr().get(),
                                  this->_output_image_adjoint_nifti_sptr->get_raw_nifti_sptr().get(),
                                  _adjoint_output_weights_sptr->get_raw_nifti_sptr().get());

    // The output should be a clone of the floating image, with data filled in from the nifti image
    output_sptr->fill(*this->_output_image_adjoint_nifti_sptr);
    this->_output_image_sptr = output_sptr;
}

namespace sirf {
template class NiftyResample<float>;
}

