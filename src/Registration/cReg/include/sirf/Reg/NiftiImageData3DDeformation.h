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
\brief Class for deformation SIRF image data.

\author Richard Brown
\author SyneRBI
*/

#pragma once

#include "sirf/Reg/NiftiImageData3DTensor.h"
#include "sirf/Reg/NonRigidTransformation.h"

namespace sirf {

// Forward declarations
template<class dataType> class NiftiImageData3D;
template<class dataType> class NiftiImageData3DDisplacement;

/*!
\ingroup Registration
\brief Class for deformation SIRF image data.

Deformation fields (as opposed to Displacement fields) describe the absolute position (in real world units) of the pixel locations on the reference image.
A deformation field of an identity transformation will contain the location of each of the pixels centroids in the world coordinates.
Here, we require ndim == 5, nt == 1 (so contains x,y,z dimensions as well as u==3 for the tensor component).
Also require intent_p1 == DEF_FIELD.

\author Richard Brown
\author SyneRBI
*/
template<class dataType>
class NiftiImageData3DDeformation : public NiftiImageData3DTensor<dataType>, public NonRigidTransformation<dataType>
{
public:
    /// Constructor
    NiftiImageData3DDeformation() {}

    /// Filename constructor
    NiftiImageData3DDeformation(const std::string &filename)
        : NiftiImageData3DTensor<dataType>(filename) { this->check_dimensions(this->_3DDef); }

    /// Nifti constructor
    NiftiImageData3DDeformation(const nifti_image &image_nifti)
        : NiftiImageData3DTensor<dataType>(image_nifti) { this->check_dimensions(this->_3DDef); }

    /// Construct from general tensor
    NiftiImageData3DDeformation(const NiftiImageData<dataType>& tensor)
        : NiftiImageData3DTensor<dataType>(tensor) { this->check_dimensions(this->_3DDef); }

    /// Construct from array
    template<class inputType>
    NiftiImageData3DDeformation(const inputType * const data, const VoxelisedGeometricalInfo3D &geom)
        : NiftiImageData3DTensor<dataType>(data, geom) { this->_nifti_image->intent_code = NIFTI_INTENT_VECTOR; this->_nifti_image->intent_p1=0; }

    /// Create from 3 individual components
    NiftiImageData3DDeformation(const NiftiImageData3D<dataType> &x, const NiftiImageData3D<dataType> &y, const NiftiImageData3D<dataType> &z)
        : NiftiImageData3DTensor<dataType>(x,y,z) { this->check_dimensions(this->_3DDef); }

    /// Create from displacement field image
    NiftiImageData3DDeformation(const NiftiImageData3DDisplacement<dataType> &disp);

    /// Create from 3D image
    void create_from_3D_image(const NiftiImageData<dataType> &image);

    /// Create from control point grid image
    void create_from_cpp(NiftiImageData3DTensor<dataType> &cpp, const NiftiImageData<dataType> &ref);

    /// Get as deformation field
    virtual NiftiImageData3DDeformation get_as_deformation_field(const NiftiImageData<dataType> &ref) const;

    /// Compose multiple transformations into single deformation field
    static NiftiImageData3DDeformation compose_single_deformation(const std::vector<const Transformation<dataType> *> &transformations, const NiftiImageData<dataType> &ref);

    /// Compose multiple transformations into single deformation field
    static NiftiImageData3DDeformation compose_single_deformation(const std::vector<std::shared_ptr<const Transformation<dataType> > > &transformations, const NiftiImageData<dataType> &ref);

    virtual ObjectHandle<DataContainer>* new_data_container_handle() const
    {
        return new ObjectHandle<DataContainer>
            (std::shared_ptr<DataContainer>(new NiftiImageData3DDeformation));
    }
    /// Write
    virtual void write(const std::string &filename) const { this->NiftiImageData<dataType>::write(filename); }
    /// Clone and return as unique pointer.
    std::unique_ptr<NiftiImageData3DDeformation> clone() const
    {
	return std::unique_ptr<NiftiImageData3DDeformation>(this->clone_impl());
    }

    /*! \brief Get inverse as unique pointer (potentially based on another image).
     *
     * Why would you want to base it on another image? Well, we might have a deformation
     * that takes us from image A to B. We'll probably want the inverse to take us from
     * image B back to A. In this case, use get_inverse(A). This is because the the deformation
     * field is defined for the reference image. In the second case, A is the reference,
     * and B is the floating image.*/
    std::unique_ptr<NiftiImageData3DDeformation> get_inverse(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr, const bool use_vtk=false) const
    {
        if (!use_vtk)
            return std::unique_ptr<NiftiImageData3DDeformation>(this->get_inverse_impl_nr(image_sptr));
        else
            return std::unique_ptr<NiftiImageData3DDeformation>(this->get_inverse_impl_vtk(image_sptr));
    }

protected:
    /// Clone helper function. Don't use.
    virtual NiftiImageData3DDeformation* clone_impl() const
    {
	return new NiftiImageData3DDeformation(*this);
    }

    /// Helper function for get_inverse (NiftyReg). Don't use.
    virtual NiftiImageData3DDeformation* get_inverse_impl_nr(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr) const;

    /// Helper function for get_inverse (VTK). Don't use.
    virtual NiftiImageData3DDeformation* get_inverse_impl_vtk(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr) const;
};
}
