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

#pragma once

#include <string>

namespace sirf {

// Forward declarations
template<class dataType> class NiftiImageData3D;
template<class dataType> class NiftiImageData3DDeformation;

/// Abstract base class for SIRFReg transformations
template<class dataType>
class Transformation
{
public:

    /// Constructor
    Transformation() {}

    /// Destructor
    virtual ~Transformation() {}

    /// Get as deformation field
    virtual NiftiImageData3DDeformation<dataType> get_as_deformation_field(const NiftiImageData3D<dataType> &ref) const = 0;

    /// Write
    virtual void write(const std::string &filename) const = 0;

protected:
    /// Check that the deformation field image matches the reference image.
    static void check_ref_and_def(const NiftiImageData3D<dataType> &ref, const NiftiImageData3DDeformation<dataType> &def);
};
}
