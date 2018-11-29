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
\brief Classes for SIRFReg transformations.
\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFREGTRANSFORMATION_H
#define _SIRFREGTRANSFORMATION_H

namespace sirf {

// Forward declarations
class NiftiImageData3D;
class NiftiImageData3DDeformation;

/// Abstract base class for SIRFReg transformations
class SIRFRegTransformation
{
public:

    /// Constructor
    SIRFRegTransformation() {}

    /// Destructor
    virtual ~SIRFRegTransformation() {}

    /// Get as deformation field
    virtual NiftiImageData3DDeformation get_as_deformation_field(const NiftiImageData3D &ref) const = 0;

protected:
    /// Check that the deformation field image matches the reference image.
    void check_ref_and_def(const NiftiImageData3D &ref, const NiftiImageData3DDeformation &def) const;
};
}

#endif
