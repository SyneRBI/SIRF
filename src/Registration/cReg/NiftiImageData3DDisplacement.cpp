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
\brief Classes for SIRFReg transformations.

\author Richard Brown
\author CCP PETMR
*/

#include "NiftiImageData3DDisplacement.h"
#include <_reg_localTrans.h>

using namespace sirf;

template<class dataType>
void NiftiImageData3DDisplacement<dataType>::create_from_def(const NiftiImageData3DDeformation<dataType> &def)
{
    // Get the disp field from the def field
    NiftiImageData3DTensor<dataType> temp = def;
    reg_getDisplacementFromDeformation(temp.get_raw_nifti_sptr().get());
    temp.get_raw_nifti_sptr()->intent_p1 = DISP_FIELD;
    *this = temp;
}

template<class dataType>
void NiftiImageData3DDisplacement<dataType>::create_from_3D_image(const NiftiImageData3D<dataType> &image)
{
    this->NiftiImageData3DTensor<dataType>::create_from_3D_image(image);
    this->_nifti_image->intent_p1 = 1;
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DDisplacement<dataType>::get_as_deformation_field(const NiftiImageData3D<dataType> &ref) const
{
    NiftiImageData3DDeformation<dataType> def;
    def.create_from_disp(*this);
    this->check_ref_and_def(ref,def);
    return def;
}

template<class dataType>
std::shared_ptr<SIRFRegTransformation<dataType> > NiftiImageData3DDisplacement<dataType>::get_clone_sptr() const
{
    return std::shared_ptr<NiftiImageData3DDisplacement<dataType> >
            (new NiftiImageData3DDisplacement<dataType>(*this));
}

namespace sirf {
template class NiftiImageData3DDisplacement<float>;
}
