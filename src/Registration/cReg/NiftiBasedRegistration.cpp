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
\brief Base class for all NIfTI-based registration algorithms.

\author Richard Brown
\author SyneRBI
*/

#include "sirf/Reg/NiftiBasedRegistration.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"

using namespace sirf;

template<class dataType>
void NiftiBasedRegistration<dataType>::convert_to_NiftiImageData_if_not_already(std::shared_ptr<const NiftiImageData3D<dataType> > &output_sptr, const std::shared_ptr<const ImageData> &input_sptr)
{
    // Try to dynamic cast from ImageData to (const) NiftiImageData. This will only succeed if original type was NiftiImageData
    output_sptr = std::dynamic_pointer_cast<const NiftiImageData3D<dataType> >(input_sptr);
    // If output is a null pointer, it means that a different image type was supplied (e.g., STIRImageData).
    // In this case, construct a NiftiImageData
    if (!output_sptr)
        output_sptr = std::make_shared<const NiftiImageData3D<dataType> >(*input_sptr);
}

template<class dataType>
const std::shared_ptr<const Transformation<dataType> > NiftiBasedRegistration<dataType>::get_displacement_field_forward_sptr(const unsigned idx) const
{
    // Get deformation as NiftiImageData3DDisplacement (from Transformation)
    std::shared_ptr<const NiftiImageData3DDeformation<dataType> > def_fwd = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<dataType> >(this->get_deformation_field_forward_sptr(idx));
    return std::move(std::make_shared<NiftiImageData3DDisplacement<dataType> >(*def_fwd));
}

template<class dataType>
const std::shared_ptr<const Transformation<dataType> > NiftiBasedRegistration<dataType>::get_displacement_field_inverse_sptr(const unsigned idx) const
{
    // Get deformation as NiftiImageData3DDisplacement (from Transformation)
    std::shared_ptr<const NiftiImageData3DDeformation<dataType> > def_inv = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<dataType> >(this->get_deformation_field_inverse_sptr(idx));
    return std::move(std::make_shared<NiftiImageData3DDisplacement<dataType> >(*def_inv));
}

namespace sirf {
template class NiftiBasedRegistration<float>;
}
