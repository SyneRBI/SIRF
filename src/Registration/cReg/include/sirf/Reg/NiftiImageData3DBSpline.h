/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 University College London

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
\brief Class for b-spline control point grid SIRF image data.

\author Richard Brown
\author SyneRBI
*/

#pragma once

#include "sirf/Reg/NiftiImageData3DTensor.h"
#include "sirf/Reg/NonRigidTransformation.h"

namespace sirf {

/*!
\ingroup Registration
\brief Class for b-spline control point grid SIRF image data.

\author Richard Brown
\author SyneRBI
*/
template<class dataType>
class NiftiImageData3DBSpline : public NiftiImageData3DTensor<dataType>, public NonRigidTransformation<dataType>
{
public:
    /// Constructor
    NiftiImageData3DBSpline() {}

    /// Filename constructor
    NiftiImageData3DBSpline(const std::string &filename)
        : NiftiImageData3DTensor<dataType>(filename) { this->check_dimensions(this->_3DBSpl); }

    /// Nifti constructor
    NiftiImageData3DBSpline(const nifti_image &image_nifti)
        : NiftiImageData3DTensor<dataType>(image_nifti) { this->check_dimensions(this->_3DBSpl); }

    /// Construct from general tensor
    NiftiImageData3DBSpline(const NiftiImageData<dataType>& tensor)
        : NiftiImageData3DTensor<dataType>(tensor) { this->check_dimensions(this->_3DBSpl); }

    /// Construct from array
    template<class inputType>
    NiftiImageData3DBSpline(const inputType * const data, const VoxelisedGeometricalInfo3D &geom)
        : NiftiImageData3DTensor<dataType>(data, geom) { this->_nifti_image->intent_code = NIFTI_INTENT_VECTOR; this->_nifti_image->intent_p1=SPLINE_VEL_GRID; }

    /// Create from 3 individual components
    NiftiImageData3DBSpline(const NiftiImageData3D<dataType> &x, const NiftiImageData3D<dataType> &y, const NiftiImageData3D<dataType> &z)
        : NiftiImageData3DTensor<dataType>(x,y,z) { this->check_dimensions(this->_3DBSpl); }

    /// Create from 3D image
    void create_from_3D_image(const NiftiImageData<dataType> &image);

    /// Get as deformation field
    virtual NiftiImageData3DDeformation<dataType> get_as_deformation_field(const NiftiImageData<dataType> &ref) const;

    /// New data handle
    virtual ObjectHandle<DataContainer>* new_data_container_handle() const
    {
        return new ObjectHandle<DataContainer>
            (std::shared_ptr<DataContainer>(new NiftiImageData3DBSpline));
    }
    /// Write
    virtual void write(const std::string &filename) const { this->NiftiImageData<dataType>::write(filename); }
    /// Clone and return as unique pointer.
    std::unique_ptr<NiftiImageData3DBSpline> clone() const
    {
        return std::unique_ptr<NiftiImageData3DBSpline>(this->clone_impl());
    }

    /*! \brief Get inverse as unique pointer (potentially based on another image).
     *
     * Why would you want to base it on another image? Well, we might have a deformation
     * that takes us from image A to B. We'll probably want the inverse to take us from
     * image B back to A. In this case, use get_inverse(A). This is because the the deformation
     * field is defined for the reference image. In the second case, A is the reference,
     * and B is the floating image.*/
    std::unique_ptr<NiftiImageData3DBSpline> get_inverse(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr, const bool use_vtk=false) const
    {
        throw std::runtime_error("NiftiImageData3DBSpline::get_inverse: not yet implemented");
    }


protected:
    /// Clone helper function. Don't use.
    virtual NiftiImageData3DBSpline* clone_impl() const
    {
        return new NiftiImageData3DBSpline(*this);
    }

    /// Helper function for get_inverse (NiftyReg). Don't use.
    virtual NiftiImageData3DBSpline* get_inverse_impl_nr(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr) const;

    /// Helper function for get_inverse (VTK). Don't use.
    virtual NiftiImageData3DBSpline* get_inverse_impl_vtk(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr) const;
};
}
