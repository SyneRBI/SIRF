/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

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
\brief Class for tensor SIRF image data.

\author Richard Brown
\author SyneRBI
*/

#pragma once

#include "sirf/Reg/NiftiImageData.h"

namespace sirf {

// Forward declarations
template<class dataType> class NiftiImageData3D;

/*!
\ingroup Registration
\brief Class for tensor SIRF image data.

Here, we require ndim == 5, nt == 1 (so contains x,y,z dimensions as well as u==3 for the tensor component).
This is the general tensor class, so we do not care if the image is a deformation or displacement field (or any other type).
Hence, any value of intent_p1 is fine.

\author Richard Brown
\author SyneRBI
*/
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

    /// Construct from array
    template<class inputType>
    NiftiImageData3DTensor(const inputType * const data, const VoxelisedGeometricalInfo3D &geom)
        : NiftiImageData<dataType>(data, geom, true) { this->_nifti_image->intent_code = NIFTI_INTENT_VECTOR; this->_nifti_image->intent_p1=-1; }

    /// Create from 3 individual components
    NiftiImageData3DTensor(const NiftiImageData3D<dataType> &x, const NiftiImageData3D<dataType> &y, const NiftiImageData3D<dataType> &z);

    /// Create from 3D image (fill with zeroes).
    virtual void create_from_3D_image(const NiftiImageData<dataType> &image);

    /// Get tensor component (x, y or z)
    std::shared_ptr<NiftiImageData<dataType> > get_tensor_component(const int component) const;

    /// Save to file as x-, y-, z-components
    void write_split_xyz_components(const std::string &filename_pattern, const int datatype=-1) const;

    /// Save to file as x-, y-, z-components
    void write_split_xyz_components(const std::string &filename_x, const std::string &filename_y, const std::string &filename_z, const int datatype=-1) const;

    /// Flip component of nu
    void flip_component(const int dim);

    /// Tensor component maths
    void tensor_component_maths(const int dim, const std::shared_ptr<const ImageData> &scalar_im_sptr, const typename NiftiImageData<dataType>::MathsType maths_type);

    /// Multiply tensor component by image
    void multiply_tensor_component(const int dim, const std::shared_ptr<const ImageData> &scalar_im_sptr);

    /// Add image to tensor component
    void add_to_tensor_component(const int dim, const std::shared_ptr<const ImageData> &scalar_im_sptr);

    virtual ObjectHandle<DataContainer>* new_data_container_handle() const
    {
        return new ObjectHandle<DataContainer>
            (std::shared_ptr<DataContainer>(new NiftiImageData3DTensor));
    }
    /// Clone and return as unique pointer.
    std::unique_ptr<NiftiImageData3DTensor> clone() const
    {
	return std::unique_ptr<NiftiImageData3DTensor>(this->clone_impl());
    }
protected:
    /// Clone helper function. Don't use.
    virtual NiftiImageData3DTensor* clone_impl() const
    {
	return new NiftiImageData3DTensor(*this);
    }
};
}
