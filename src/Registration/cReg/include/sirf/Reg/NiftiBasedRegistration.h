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
\brief Base class for all NIfTI-based registrations.

\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include "sirf/Reg/Registration.h"

namespace sirf {

/// Forward declarations
template<class dataType> class NiftiImageData3D;

/*!
\ingroup Registration
\brief Base class for all NIfTI-based registrations.

\author Richard Brown
\author CCP PETMR
*/
template<class dataType>
class NiftiBasedRegistration : public Registration<dataType>
{
public:

    /// Constructor
    NiftiBasedRegistration() {}

    /// Destructor
    virtual ~NiftiBasedRegistration() {}

    /// Get forward deformation field image
    virtual const std::shared_ptr<const Transformation<dataType> > get_deformation_field_forward_sptr(const unsigned idx = 0) const;

    /// Get inverse deformation field image
    virtual const std::shared_ptr<const Transformation<dataType> > get_deformation_field_inverse_sptr(const unsigned idx = 0) const;

    /// Get registered image as NiftiImageData3D
    const std::shared_ptr<const NiftiImageData3D<dataType> > get_output_sptr(const unsigned idx = 0) const { return _warped_images_nifti.at(idx); }

    /// Convert an ImageData to NiftiImageData. Try to dynamic pointer cast, else create new image.
    static void convert_to_NiftiImageData_if_not_already(std::shared_ptr<const NiftiImageData3D<dataType> > &output_sptr, const std::shared_ptr<const ImageData> &input_sptr);

protected:

    /// Reference image (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _reference_image_nifti_sptr;
    /// Floating image (as NiftiImageData3D)
    std::vector<std::shared_ptr<const NiftiImageData3D<dataType> > > _floating_images_nifti;
    /// Output (as NiftiImageData3D)
    std::vector<std::shared_ptr<NiftiImageData3D<dataType> > > _warped_images_nifti;
};
}
