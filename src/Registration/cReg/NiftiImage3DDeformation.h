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

#ifndef _NIFTIIMAGE3DDEFORMATION_H_
#define _NIFTIIMAGE3DDEFORMATION_H_

#include "NiftiImage3DTensor.h"
#include <_reg_maths.h>

namespace sirf {
class NiftiImage3D;

/// SIRF nifti image data deformation field image
class NiftiImage3DDeformation : public NiftiImage3DTensor
{
public:
    /// Constructor
    NiftiImage3DDeformation() {}

    /// Filename constructor
    NiftiImage3DDeformation(const std::string &filename)
        : NiftiImage3DTensor(filename) { check_dimensions(_3DDef); }

    /// Nifti constructor
    NiftiImage3DDeformation(const nifti_image &image_nifti)
        : NiftiImage3DTensor(image_nifti) { check_dimensions(_3DDef); }

    /// Nifti shared_ptr constructor
    NiftiImage3DDeformation(const std::shared_ptr<nifti_image> image_nifti)
        : NiftiImage3DTensor(image_nifti) { check_dimensions(_3DDef); }

    /// Construct from general tensor
    NiftiImage3DDeformation(const NiftiImage& tensor)
        : NiftiImage3DTensor(tensor) { check_dimensions(_3DDef); }

    /// Deep copy
    NiftiImage3DDeformation deep_copy() const
    { return this->NiftiImage3DTensor::deep_copy(); }

    /// Create from 3D image
    void create_from_3D_image(const NiftiImage3D &image)
    {
        this->NiftiImage3DTensor::create_from_3D_image(image);
        _nifti_image->intent_p1 = 0;
    }
};
}

#endif
