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
\brief Class for deformation/displacement SIRF image data.

\author Richard Brown
\author CCP PETMR
*/

#ifndef _NIFTIIMAGE3DDISPLACEMENT_H_
#define _NIFTIIMAGE3DDISPLACEMENT_H_

#include "NiftiImage3DTensor.h"
#include "NiftiImage3DDeformation.h"
#include "SIRFRegTransformation.h"
#include <_reg_maths.h>

namespace sirf {
class NiftiImage3D;

/// SIRF nifti image data displacement field image
class NiftiImage3DDisplacement : public NiftiImage3DTensor, public SIRFRegTransformation
{
public:
    /// Constructor
    NiftiImage3DDisplacement() {}

    /// Filename constructor
    NiftiImage3DDisplacement(const std::string &filename)
        : NiftiImage3DTensor(filename) { check_dimensions(_3DDisp); }

    /// Nifti constructor
    NiftiImage3DDisplacement(const nifti_image &image_nifti)
        : NiftiImage3DTensor(image_nifti) { check_dimensions(_3DDisp); }

    /// Nifti shared_ptr constructor
    NiftiImage3DDisplacement(const std::shared_ptr<nifti_image> image_nifti)
        : NiftiImage3DTensor(image_nifti) { check_dimensions(_3DDisp); }

    /// Construct from general tensor
    NiftiImage3DDisplacement(const NiftiImage& tensor)
        : NiftiImage3DTensor(tensor) { check_dimensions(_3DDisp); }

    /// Create from 3 individual components
    NiftiImage3DDisplacement(const NiftiImage3D &x, const NiftiImage3D &y, const NiftiImage3D &z)
        : NiftiImage3DTensor(x,y,z) { _nifti_image->intent_p1 = 1; }

    /// Create from deformation field image
    void create_from_def(const NiftiImage3DDeformation &im);

    /// Deep copy
    NiftiImage3DDisplacement deep_copy() const
    { return this->NiftiImage::deep_copy(); }

    /// Create from 3D image
    void create_from_3D_image(const NiftiImage3D &image);

    /// Get as deformation field
    virtual NiftiImage3DDeformation get_as_deformation_field(const NiftiImage3D &ref) const;
};
}

#endif
