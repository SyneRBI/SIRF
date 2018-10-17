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

#ifndef _NIFTIIMAGEDATA3DDISPLACEMENT_H_
#define _NIFTIIMAGEDATA3DDISPLACEMENT_H_

#include "NiftiImageData3DTensor.h"
#include "NiftiImageData3DDeformation.h"
#include "SIRFRegTransformation.h"
#include <_reg_maths.h>

namespace sirf {
class NiftiImageData3D;

/// SIRF nifti image data displacement field image
class NiftiImageData3DDisplacement : public NiftiImageData3DTensor, public SIRFRegTransformation
{
public:
    /// Constructor
    NiftiImageData3DDisplacement() {}

    /// Filename constructor
    NiftiImageData3DDisplacement(const std::string &filename)
        : NiftiImageData3DTensor(filename) { check_dimensions(_3DDisp); }

    /// Nifti constructor
    NiftiImageData3DDisplacement(const nifti_image &image_nifti)
        : NiftiImageData3DTensor(image_nifti) { check_dimensions(_3DDisp); }

    /// Nifti std::shared_ptr constructor
    NiftiImageData3DDisplacement(const std::shared_ptr<nifti_image> image_nifti)
        : NiftiImageData3DTensor(image_nifti) { check_dimensions(_3DDisp); }

    /// Construct from general tensor
    NiftiImageData3DDisplacement(const NiftiImageData& tensor)
        : NiftiImageData3DTensor(tensor) { check_dimensions(_3DDisp); }

    /// Create from 3 individual components
    NiftiImageData3DDisplacement(const NiftiImageData3D &x, const NiftiImageData3D &y, const NiftiImageData3D &z)
        : NiftiImageData3DTensor(x,y,z) { _nifti_image->intent_p1 = 1; }

    /// Create from deformation field image
    void create_from_def(const NiftiImageData3DDeformation &im);

    /// Deep copy
    NiftiImageData3DDisplacement deep_copy() const
    { return this->NiftiImageData::deep_copy(); }

    /// Create from 3D image
    void create_from_3D_image(const NiftiImageData3D &image);

    /// Get as deformation field
    virtual NiftiImageData3DDeformation get_as_deformation_field(const NiftiImageData3D &ref) const;
};
}

#endif
