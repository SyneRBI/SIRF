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

#ifndef _NIFTIIMAGEDATA3DDEFORMATION_H_
#define _NIFTIIMAGEDATA3DDEFORMATION_H_

#include "NiftiImageData3DTensor.h"
#include "SIRFRegTransformation.h"
#include <_reg_maths.h>

namespace sirf {
class NiftiImageData3D;

/// SIRF nifti image data deformation field image
class NiftiImageData3DDeformation : public NiftiImageData3DTensor, public SIRFRegTransformation
{
public:
    /// Constructor
    NiftiImageData3DDeformation() {}

    /// Filename constructor
    NiftiImageData3DDeformation(const std::string &filename)
        : NiftiImageData3DTensor(filename) { check_dimensions(_3DDef); }

    /// Nifti constructor
    NiftiImageData3DDeformation(const nifti_image &image_nifti)
        : NiftiImageData3DTensor(image_nifti) { check_dimensions(_3DDef); }

    /// Nifti std::shared_ptr constructor
    NiftiImageData3DDeformation(const std::shared_ptr<nifti_image> image_nifti)
        : NiftiImageData3DTensor(image_nifti) { check_dimensions(_3DDef); }

    /// Construct from general tensor
    NiftiImageData3DDeformation(const NiftiImageData& tensor)
        : NiftiImageData3DTensor(tensor) { check_dimensions(_3DDef); }

    /// Create from 3 individual components
    NiftiImageData3DDeformation(const NiftiImageData3D &x, const NiftiImageData3D &y, const NiftiImageData3D &z)
        : NiftiImageData3DTensor(x,y,z) { check_dimensions(_3DDef); }

    /// Create from deformation field image
    void create_from_disp(const NiftiImageData3DDisplacement &im);

    /// Deep copy
    NiftiImageData3DDeformation deep_copy() const
    { return this->NiftiImageData3DTensor::deep_copy(); }

    /// Create from 3D image
    void create_from_3D_image(const NiftiImageData3D &image);

    /// Create from CPP image
    void create_from_cpp(NiftiImageData3DTensor &cpp, const NiftiImageData3D &ref);

    /// Get as deformation field
    virtual NiftiImageData3DDeformation get_as_deformation_field(const NiftiImageData3D &ref) const;

    /// Compose multiple transformations into single deformation field
    static NiftiImageData3DDeformation compose_single_deformation(const std::vector<SIRFRegTransformation*> &transformations, const NiftiImageData3D &ref);

    /// Compose multiple transformations into single deformation field
    static NiftiImageData3DDeformation compose_single_deformation(const std::vector<std::shared_ptr<SIRFRegTransformation> > &transformations, const NiftiImageData3D &ref);
};
}

#endif
