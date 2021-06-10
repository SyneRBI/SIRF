/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2019 University College London

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
\brief Base class for transformations.
\author Richard Brown
\author SyneRBI
*/

#pragma once

#include <string>

namespace sirf {

// Forward declarations
template<class dataType> class NiftiImageData;
template<class dataType> class NiftiImageData3D;
template<class dataType> class NiftiImageData3DDeformation;

/*!
\ingroup Registration
\brief Base class for transformations.

All transformations need to be able to convert themselves
to a deformation field. In this fashion, they can be composed
into a single transformation.

\author Richard Brown
\author SyneRBI
*/template<class dataType>
class Transformation
{
public:

    /// Constructor
    Transformation() {}

    /// Destructor
    virtual ~Transformation() {}

    /// Get as deformation field.
    ///
    /// Reference image **must** be used when converting a transformation matrix to a deformation.
    /// For displacements and deformations, the reference can be used optionally. It **should** be used when composing transformations to be used for resampling
    /// But is probably unnecessary for simply concatenating deformations.
    virtual NiftiImageData3DDeformation<dataType> get_as_deformation_field(const NiftiImageData<dataType> &ref, const bool use_ref = true) const = 0;

    /// Write
    virtual void write(const std::string &filename) const = 0;

protected:
    /// Check that the deformation field image matches the reference image.
    static void check_ref_and_def(const NiftiImageData<dataType> &ref, const NiftiImageData3DDeformation<dataType> &def);
};
}
