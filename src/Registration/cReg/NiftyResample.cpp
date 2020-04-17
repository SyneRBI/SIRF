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
using namespace detail;

template<class dataType>
static void convert_ImageData_to_ComplexNiftiImageData(ComplexNiftiImageData<dataType> &output, const std::shared_ptr<const ImageData> input_sptr)
{
    // if input is real, only convert first bit
    if (!input_sptr->is_complex()) {
        output.real() = std::make_shared<NiftiImageData<dataType> >(*input_sptr);
        output.imag().reset();
    }
    // if input is complex, only set both
    else {
        std::shared_ptr<NiftiImageData<dataType> > &output_real = output.real();
        std::shared_ptr<NiftiImageData<dataType> > &output_imag = output.imag();
        NiftiImageData<dataType>::construct_NiftiImageData_from_complex_im(output_real,output_imag,input_sptr);
    }
}

template<class dataType>
static void set_up_output_image(ComplexNiftiImageData<dataType> &output,
                                const ComplexNiftiImageData<dataType> im_for_shape,
                                const ComplexNiftiImageData<dataType> im_for_metadata)
{
    // The output is a mixture between the reference and floating images.
    const nifti_image * const im_for_shape_ptr = im_for_shape.real()->get_raw_nifti_sptr().get();
    const nifti_image * const im_for_metadata_ptr = im_for_metadata.real()->get_raw_nifti_sptr().get();

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
    // If the im_for_shape is complex (i.e., has two components), then output should be the same
    output.real() = std::make_shared<NiftiImageData<dataType> >(*output_ptr);
    if (im_for_shape.is_complex())
        output.imag() = std::make_shared<NiftiImageData<dataType> >(*output_ptr);

    // Delete the original pointer
    nifti_image_free(output_ptr);
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
    // Use the reference regardless of forward/adjoint.
    this->_deformation_sptr = std::make_shared<NiftiImageData3DDeformation<dataType> >(
            NiftiImageData3DDeformation<dataType>::compose_single_deformation(
                this->_transformations,*this->_reference_image_niftis.real()));

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
    set_up_output_image(_output_image_forward_niftis, _reference_image_niftis, _floating_image_niftis);

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
    set_up_output_image(_output_image_adjoint_niftis, _floating_image_niftis, _reference_image_niftis);

    // Get deformation spacing
    float control_point_grid_spacing[3] = {
        _deformation_sptr->get_raw_nifti_sptr()->dx,
        _deformation_sptr->get_raw_nifti_sptr()->dy,
        _deformation_sptr->get_raw_nifti_sptr()->dz
    };

    // Get the deformation field
    nifti_image *def_ptr = _deformation_sptr->get_raw_nifti_sptr().get();
    // Copy of reference image
    NiftiImageData<dataType> ref = *this->_reference_image_niftis.real();
    nifti_image *ref_ptr = ref.get_raw_nifti_sptr().get();

    _adjoint_transformer_sptr = std::make_shared<NiftyMoMo::BSplineTransformation>(
                ref_ptr,
                /*num levels to perform*/1U,
                control_point_grid_spacing);

    _adjoint_transformer_sptr->set_interpolation(this->_interpolation_type);
    _adjoint_transformer_sptr->SetPaddingValue(this->_padding_value);
    _adjoint_transformer_sptr->setDVF(def_ptr);

    // Set up the input and output weights
    _adjoint_input_weights_sptr = this->_reference_image_niftis.real()->clone();
    _adjoint_input_weights_sptr->fill(1.f);
    _adjoint_output_weights_sptr = this->_output_image_adjoint_niftis.real()->clone();

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
    convert_ImageData_to_ComplexNiftiImageData(this->_reference_image_niftis, this->_reference_image_sptr);
    convert_ImageData_to_ComplexNiftiImageData(this->_floating_image_niftis,  this->_floating_image_sptr);

#ifndef NDEBUG
    if (!_reference_image_niftis.is_complex() && !_floating_image_niftis.is_complex())
        std::cout << "\nNiftyResample: Reference = real, floating = real. Output forward will = real, output adjoint will = real.\n";
    else if (_reference_image_niftis.is_complex() && _floating_image_niftis.is_complex())
        std::cout << "\nNiftyResample: Reference = complex, floating = complex. Output forward will = complex, output adjoint will = complex.\n";
    else if (_reference_image_niftis.is_complex())
        std::cout << "\nNiftyResample: Reference = complex, floating = real. Output forward will = complex (with 0 imaginary), output adjoint will = real (only real component of reference will be resampled).\n";
    else if (_floating_image_niftis.is_complex())
        std::cout << "\nNiftyResample: Reference = real, floating = complex. Output forward will = real (only real component of floating will be resampled), output adjoint will = complex (with 0 imaginary).\n";
    else
        throw std::runtime_error("Shouldn't be here");
#endif
}

template<class dataType>
static void check_images_match(
        const ComplexNiftiImageData<dataType> im1,
        const ComplexNiftiImageData<dataType> im2,
        const std::string &explanation)
{
    if (im1.is_complex() != im2.is_complex())
        throw std::runtime_error(explanation);
    for (unsigned i=0; i<im1.size(); ++i)
        if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(*im1.at(i),*im2.at(i),false))
            throw std::runtime_error(explanation);
}

template<class dataType>
static void set_post_resample_outputs(std::shared_ptr<ImageData> &output_to_return_sptr, std::shared_ptr<ImageData> &output_as_member_sptr, const ComplexNiftiImageData<dataType> resampled_niftis)
{
    // If output is only real, set that
    if (!output_to_return_sptr->is_complex())
        output_to_return_sptr->fill(*resampled_niftis.real());
    // Else, set the complex bit
    else {
        NumberType::Type output_num_type = (*output_to_return_sptr->begin()).get_typeID();
        if (output_num_type != NumberType::CXFLOAT)
            throw std::runtime_error("NiftyResample: Only complex type currently supported is complex float");
        ImageData::Iterator &it_out = output_to_return_sptr->begin();
        auto &it_real = resampled_niftis.real()->begin();
        auto &it_imag = resampled_niftis.imag()->begin();
        for (; it_out!=output_to_return_sptr->end(); ++it_real, ++it_imag, ++it_out) {
            complex_float_t cmplx_flt(*it_real,*it_imag);
            *it_out = NumRef((void *)&cmplx_flt, output_num_type);
        }
    }

    // Copy the output so that backwards compatibility of get_output() is preserved.
    output_as_member_sptr = output_to_return_sptr;
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
    ComplexNiftiImageData<dataType> input_niftis, output_niftis;
    convert_ImageData_to_ComplexNiftiImageData(input_niftis, input_sptr);
    convert_ImageData_to_ComplexNiftiImageData(output_niftis, output_sptr);

    // Check that the metadata match
    check_images_match(input_niftis, this->_floating_image_niftis,
                       "NiftyResample::forward: Metadata of input image should match floating image.");
    check_images_match(output_niftis, this->_reference_image_niftis,
                       "NiftyResample::forward: Metadata of output image should match reference image.");

    // Loop over as many output images
    for (unsigned i=0; i<_output_image_forward_niftis.size(); ++i) {

        reg_resampleImage(input_niftis.at(i)->get_raw_nifti_sptr().get(),
                          _output_image_forward_niftis.at(i)->get_raw_nifti_sptr().get(),
                          _deformation_sptr->get_raw_nifti_sptr().get(),
                          nullptr,
                          this->_interpolation_type,
                          this->_padding_value);
    }

    set_post_resample_outputs(output_sptr, this->_output_image_sptr, _output_image_forward_niftis);
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
    ComplexNiftiImageData<dataType> input_niftis, output_niftis;
    convert_ImageData_to_ComplexNiftiImageData(input_niftis, input_sptr);
    convert_ImageData_to_ComplexNiftiImageData(output_niftis, output_sptr);

    // Check that the metadata match
    check_images_match(input_niftis, this->_reference_image_niftis,
                       "NiftyResample::adjoint: Metadata of input image should match reference image.");
    check_images_match(output_niftis, this->_floating_image_niftis,
                       "NiftyResample::adjoint: Metadata of output image should match floating image.");

    // Loop over the real and potentially imaginary parts
    for (unsigned i=0; i<output_niftis.size(); ++i) {

        // set output to zero
        this->_output_image_adjoint_niftis.at(i)->fill(0.f);

        _adjoint_transformer_sptr->
                TransformImageAdjoint(input_niftis.at(i)->get_raw_nifti_sptr().get(),
                                      _adjoint_input_weights_sptr->get_raw_nifti_sptr().get(),
                                      this->_output_image_adjoint_niftis.at(i)->get_raw_nifti_sptr().get(),
                                      _adjoint_output_weights_sptr->get_raw_nifti_sptr().get());
    }

    set_post_resample_outputs(output_sptr, this->_output_image_sptr, _output_image_adjoint_niftis);
}

namespace sirf {
template class NiftyResample<float>;
}

