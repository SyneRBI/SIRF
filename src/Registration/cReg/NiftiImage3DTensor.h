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

#ifndef _NIFTIIMAGE3DTENSOR_H_
#define _NIFTIIMAGE3DTENSOR_H_

#include "NiftiImage.h"
#include <_reg_maths.h>

namespace sirf {
class NiftiImage3D;

/// SIRF image data
class NiftiImage3DTensor : public NiftiImage
{
public:

    /// Constructor
    NiftiImage3DTensor() {}

    /// Construct 3D from general case
    NiftiImage3DTensor(const NiftiImage& general)
        : NiftiImage(general) { check_dimensions(_3DTensor); }

    /// Filename constructor
    NiftiImage3DTensor(const std::string &filename)
        : NiftiImage(filename) { check_dimensions(_3DTensor); }

    /// Nifti constructor
    NiftiImage3DTensor(const nifti_image &image_nifti)
        : NiftiImage(image_nifti) { check_dimensions(_3DTensor); }

    /// Nifti shared_ptr constructor
    NiftiImage3DTensor(const std::shared_ptr<nifti_image> image_nifti)
        : NiftiImage(image_nifti) { check_dimensions(_3DTensor); }

    /// Create from 3 individual components
    NiftiImage3DTensor(const NiftiImage3D &x, const NiftiImage3D &y, const NiftiImage3D &z);

    /// Create from 3D image.
    virtual void create_from_3D_image(const NiftiImage3D &image);

    /// Save to file as x-, y-, z-components
    void save_to_file_split_xyz_components(const std::string &filename_pattern, const int datatype=-1) const;

    /// Save to file as x-, y-, z-components
    void save_to_file_split_xyz_components(const std::string &filename_x, const std::string &filename_y, const std::string &filename_z, const int datatype=-1) const;

    /// Flip component of nu
    void flip_component(const int dim);

    /// Deep copy
    NiftiImage3DTensor deep_copy() const
    { return this->NiftiImage::deep_copy(); }
};
}

#endif
