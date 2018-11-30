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

#ifndef _NIFTIIMAGEDATA3DTENSOR_H_
#define _NIFTIIMAGEDATA3DTENSOR_H_

#include "NiftiImageData.h"
#include <_reg_maths.h>

namespace sirf {

// Forward declarations
template<class dataType> class NiftiImageData3D;

/// SIRF image data
template<class dataType>
class NiftiImageData3DTensor : public NiftiImageData<dataType>
{
public:

    /// Constructor
    NiftiImageData3DTensor() {}

    /// Construct 3D from general case
    NiftiImageData3DTensor(const NiftiImageData<dataType>& general)
        : NiftiImageData<dataType>(general) { this->check_dimensions(this->_3DTensor); }

    /// Filename constructor
    NiftiImageData3DTensor(const std::string &filename)
        : NiftiImageData<dataType>(filename) { this->check_dimensions(this->_3DTensor); }

    /// Nifti constructor
    NiftiImageData3DTensor(const nifti_image &image_nifti)
        : NiftiImageData<dataType>(image_nifti) { this->check_dimensions(this->_3DTensor); }

    /// Nifti std::shared_ptr constructor
    NiftiImageData3DTensor(const std::shared_ptr<nifti_image> image_nifti)
        : NiftiImageData<dataType>(image_nifti) { this->check_dimensions(this->_3DTensor); }

    /// Create from 3 individual components
    NiftiImageData3DTensor(const NiftiImageData3D<dataType> &x, const NiftiImageData3D<dataType> &y, const NiftiImageData3D<dataType> &z);

    /// Create from 3D image.
    virtual void create_from_3D_image(const NiftiImageData3D<dataType> &image);

    /// Save to file as x-, y-, z-components
    void save_to_file_split_xyz_components(const std::string &filename_pattern, const int datatype=-1) const;

    /// Save to file as x-, y-, z-components
    void save_to_file_split_xyz_components(const std::string &filename_x, const std::string &filename_y, const std::string &filename_z, const int datatype=-1) const;

    /// Flip component of nu
    void flip_component(const int dim);

    /// Deep copy
    NiftiImageData3DTensor deep_copy() const
    { return this->NiftiImageData<dataType>::deep_copy(); }
};
}

#endif
