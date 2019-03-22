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
\brief Base class for all NiftyReg registrations.

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
\brief Base class for all NiftyReg registrations.

\author Richard Brown
\author CCP PETMR
*/
template<class dataType>
class NiftyRegistration : public Registration<dataType>
{
public:

    /// Constructor
    NiftyRegistration() {}

    /// Destructor
    virtual ~NiftyRegistration() {}

    /// Get forward deformation field image
    virtual const std::shared_ptr<const Transformation<dataType> > get_deformation_field_forward_sptr() const;

    /// Get inverse deformation field image
    virtual const std::shared_ptr<const Transformation<dataType> > get_deformation_field_inverse_sptr() const;

    /// Get registered image as NiftiImageData3D
    const std::shared_ptr<const NiftiImageData3D<dataType> > get_output_sptr() const { return _warped_image_nifti_sptr; }

protected:

    /// Set up inputs
    void set_up_inputs();

    /// Reference image (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _reference_image_nifti_sptr;
    /// Floating image (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _floating_image_nifti_sptr;
    /// Floating mask (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _floating_mask_nifti_sptr;
    /// Reference mask (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _reference_mask_nifti_sptr;
    /// Output (as NiftiImageData3D)
    std::shared_ptr<NiftiImageData3D<dataType> >       _warped_image_nifti_sptr;
};
}
