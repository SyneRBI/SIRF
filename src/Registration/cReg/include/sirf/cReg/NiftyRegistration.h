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

The parameter file should have interfile-like syntax.
The variables will be stored as a vector of floats and converted into the required type (int, unsigned int, etc) if necessary.
Multiple variables for a given parameter should be comma separated.
Spaces and tabs will be ignored.
For the title, it doesn't matter what is written as it will be ignored, but something has to be there (otherwise the first parameter will be ignored).
Possible parameters are all the Set<something> methods for each class (e.g., nifty_aladin::SetPerformRigid) and should be written in the parameter file without the "Set" (e.g., PerformRigid).

An example is given below:
    SomeTitle :=
        ReferenceTimePoint := 1
        FloatingTimePoint := 2
        LinearEnergyWeights := 1.5,1
        AdditiveMC :=
    end :=

More examples can be found in data/examples/Registration/paramFiles

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFREGNIFTYREGISTRATION_H_
#define _SIRFREGNIFTYREGISTRATION_H_

#include "sirf/cReg/Registration.h"

namespace sirf {

/// Forward declarations
template<class dataType> class NiftiImageData3D;

/// Base class for registration algorithms wrapped by SIRFReg
template<class dataType>
class SIRFRegNiftyRegistration : public SIRFReg<dataType>
{
public:

    /// Constructor
    SIRFRegNiftyRegistration() {}

    /// Destructor
    virtual ~SIRFRegNiftyRegistration() {}

    /// Get forward deformation field image
    virtual const std::shared_ptr<const SIRFRegTransformation<dataType> > get_deformation_field_forward() const;

    /// Get inverse deformation field image
    virtual const std::shared_ptr<const SIRFRegTransformation<dataType> > get_deformation_field_inverse() const;

    /// Get registered image as NiftiImageData3D
    const std::shared_ptr<const NiftiImageData3D<dataType> > get_output() const { return _warped_image_nifti_sptr; }

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

#endif
