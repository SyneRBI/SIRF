/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

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
\brief Base class for all SIRF registration.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/cReg/SIRFRegNiftyRegistration.h"
#include "sirf/cReg/NiftiImageData3D.h"
#include "sirf/cReg/NiftiImageData3DDisplacement.h"
#include "sirf/cReg/NiftiImageData3DDeformation.h"

using namespace sirf;

template<class dataType>
const std::shared_ptr<const SIRFRegTransformation<dataType> > SIRFRegNiftyRegistration<dataType>::get_deformation_field_forward() const
{
    // Get displacement as NiftiImageData3DDisplacement (from SIRFRegTransformation)
    std::shared_ptr<const NiftiImageData3DDisplacement<dataType> > disp_fwd = std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<dataType> >(this->_disp_image_forward_sptr);
    std::shared_ptr<NiftiImageData3DDeformation<dataType> > def_fwd = std::make_shared<NiftiImageData3DDeformation<dataType> >();
    def_fwd->create_from_disp(*disp_fwd);
    return def_fwd;
}

template<class dataType>
const std::shared_ptr<const SIRFRegTransformation<dataType> > SIRFRegNiftyRegistration<dataType>::get_deformation_field_inverse() const
{
    // Get displacement as NiftiImageData3DDisplacement (from SIRFRegTransformation)
    std::shared_ptr<const NiftiImageData3DDisplacement<dataType> > disp_inv= std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<dataType> >(this->_disp_image_inverse_sptr);
    std::shared_ptr<NiftiImageData3DDeformation<dataType> > def_inv = std::make_shared<NiftiImageData3DDeformation<dataType> >();
    def_inv->create_from_disp(*disp_inv);
    return def_inv;
}

template<class dataType>
void SIRFRegNiftyRegistration<dataType>::set_up_inputs()
{
    // Try to dynamic cast from ImageData to NiftiImageData3D. This will only succeed if original type was NiftiImageData3D
    // If the result is a null pointer, it means that a different image type was supplied (e.g., STIRImageData).
    // In this case, construct a NiftiImageData3D

    // Reference image
    this->_reference_image_nifti_sptr = std::dynamic_pointer_cast<const NiftiImageData3D<dataType> >(this->_reference_image_sptr);
    if (!this->_reference_image_nifti_sptr)
        this->_reference_image_nifti_sptr = std::make_shared<const NiftiImageData3D<dataType> >(*this->_reference_image_sptr);

    // Floating image
    this->_floating_image_nifti_sptr  = std::dynamic_pointer_cast<const NiftiImageData3D<dataType> >(this->_floating_image_sptr);
    if (!this->_floating_image_nifti_sptr)
        this->_floating_image_nifti_sptr = std::make_shared<const NiftiImageData3D<dataType> >(*this->_floating_image_sptr);

    // Reference mask (if supplied)
    if (this->_reference_mask_sptr) {
        this->_reference_mask_nifti_sptr = std::dynamic_pointer_cast<const NiftiImageData3D<dataType> >(this->_reference_mask_sptr);
        if (!this->_reference_mask_nifti_sptr)
            this->_reference_mask_nifti_sptr = std::make_shared<const NiftiImageData3D<dataType> >(*this->_reference_mask_sptr);
    }

    // Floating mask (if supplied)
    if (this->_floating_mask_sptr) {
        this->_floating_mask_nifti_sptr = std::dynamic_pointer_cast<const NiftiImageData3D<dataType> >(this->_floating_mask_sptr);
        if (!this->_floating_mask_nifti_sptr)
            this->_floating_mask_nifti_sptr = std::make_shared<const NiftiImageData3D<dataType> >(*this->_floating_mask_sptr);
    }
}

namespace sirf {
template class SIRFRegNiftyRegistration<float>;
}
