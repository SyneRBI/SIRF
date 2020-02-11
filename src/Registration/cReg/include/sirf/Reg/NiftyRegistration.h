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

#include "sirf/Reg/NiftiBasedRegistration.h"

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
class NiftyRegistration : public NiftiBasedRegistration<dataType>
{
public:

    /// Constructor
    NiftyRegistration() {}

    /// Destructor
    virtual ~NiftyRegistration() {}

protected:

    /// Set up inputs
    void set_up_inputs();

    /// Floating mask (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _floating_mask_nifti_sptr;
    /// Reference mask (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _reference_mask_nifti_sptr;
};
}
