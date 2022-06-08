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
\brief Base class for transformations.
\author Richard Brown
\author SyneRBI
*/

#pragma once

#include "sirf/Reg/Transformation.h"

namespace sirf {

/*!
\ingroup Registration
\brief Base class for non-rigid transformations.

\author Richard Brown
\author SyneRBI
*/template<class dataType>
class NonRigidTransformation : public Transformation<dataType>
{
public:

    /// Destructor
    virtual ~NonRigidTransformation() {}

#ifndef _MSC_VER // need to disable on Visual Studio for https://github.com/SyneRBI/SIRF/issues/665

    /*! \brief Get inverse as unique pointer (potentially based on another image).
     *
     * Why would you want to base it on another image? Well, we might have a deformation
     * that takes us from image A to B. We'll probably want the inverse to take us from
     * image B back to A. In this case, use get_inverse(A). This is because the the deformation
     * field is defined for the reference image. In the second case, A is the reference,
     * and B is the floating image.*/
    std::unique_ptr<NonRigidTransformation> get_inverse(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr, const bool use_vtk=false) const
    {
        if (!use_vtk)
            return std::unique_ptr<NonRigidTransformation>(this->get_inverse_impl_nr(image_sptr));
        else
            return std::unique_ptr<NonRigidTransformation>(this->get_inverse_impl_vtk(image_sptr));
    }

protected:
    /// Helper function for get_inverse (NiftyReg). Don't use.
    virtual NonRigidTransformation* get_inverse_impl_nr(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr) const = 0;

    /// Helper function for get_inverse (VTK). Don't use.
    virtual NonRigidTransformation* get_inverse_impl_vtk(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr = nullptr) const = 0;
#endif // _MSC_VER
};

}
