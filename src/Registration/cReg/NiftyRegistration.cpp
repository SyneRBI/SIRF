/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

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
\brief Base class for all NiftyReg registration algorithms.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftyRegistration.h"

using namespace sirf;

template<class dataType>
NiftyRegistration<dataType>::NiftyRegistration()
{
    this->_warped_images.resize(1);
    this->_warped_images_nifti.resize(1);
    this->_disp_fwd_images.resize(1);
    this->_disp_inv_images.resize(1);
}

template<class dataType>
void NiftyRegistration<dataType>::set_parameter(const std::string &par, const std::string &arg1, const std::string &arg2)
{
    _extra_params.push_back(par);
    _extra_params.push_back(arg1);
    _extra_params.push_back(arg2);
}

template<class dataType>
void NiftyRegistration<dataType>::set_up_inputs()
{
    if (this->_floating_images.size()!=1)
        throw std::runtime_error("NiftyReg only accepts one floating image.");
    this->_floating_images_nifti.resize(1);

    // Try to dynamic cast from ImageData to NiftiImageData3D. This will only succeed if original type was NiftiImageData3D
    // If the result is a null pointer, it means that a different image type was supplied (e.g., STIRImageData).
    // In this case, construct a NiftiImageData3D

    // Reference and floating images
    NiftiBasedRegistration<dataType>::convert_to_NiftiImageData_if_not_already(this->_reference_image_nifti_sptr, this->_reference_image_sptr);
    NiftiBasedRegistration<dataType>::convert_to_NiftiImageData_if_not_already(this->_floating_images_nifti.at(0), this->_floating_images.at(0));

    // Reference and floating masks (if supplied)
    if (this->_reference_mask_sptr)
        NiftiBasedRegistration<dataType>::convert_to_NiftiImageData_if_not_already(this->_reference_mask_nifti_sptr, this->_reference_mask_sptr);
    if (this->_floating_mask_sptr)
        NiftiBasedRegistration<dataType>::convert_to_NiftiImageData_if_not_already(this->_floating_mask_nifti_sptr, this->_floating_mask_sptr);
}

namespace sirf {
template class NiftyRegistration<float>;
}
