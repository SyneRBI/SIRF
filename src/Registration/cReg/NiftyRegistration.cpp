/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

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
\ingroup Registration
\brief Base class for all NiftyReg registration algorithms.

\author Richard Brown
\author SyneRBI
*/

#include "sirf/Reg/NiftyRegistration.h"
#include "sirf/Reg/NiftiImageData3D.h"

using namespace sirf;

template<class dataType>
NiftyRegistration<dataType>::NiftyRegistration()
{
    this->_warped_images.resize(1);
    this->_warped_images_nifti.resize(1);
    this->_def_fwd_images.resize(1);
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
    if (this->_floating_images.size()+this->_floating_image_filenames.size() != 1)
        throw std::runtime_error("NiftyReg only accepts one floating image.");
    this->_floating_images_nifti.resize(1);

    // For reference and floating image.
    // If filename has been set, read the image.
    // But if it's been set as an ImageData, convert it

    // If image has been read via filename, read it.
    if (!this->_reference_image_filename.empty())
        this->_reference_image_nifti_sptr = std::make_shared<const NiftiImageData3D<dataType> >(this->_reference_image_filename);
    // Else, convert it
    else
        NiftiBasedRegistration<dataType>::convert_to_NiftiImageData_if_not_already(this->_reference_image_nifti_sptr, this->_reference_image_sptr);

    // If image has been read via filename, read it.
    if (this->_floating_image_filenames.size() == 1)
        this->_floating_images_nifti.at(0) = std::make_shared<const NiftiImageData3D<dataType> >(this->_floating_image_filenames.at(0));
    // Else, convert it
    else
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
