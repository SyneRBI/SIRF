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

#include "NiftiImage3DDisplacement.h"
#include <_reg_localTrans.h>

using namespace std;
using namespace sirf;

void NiftiImage3DDisplacement::create_from_def(const NiftiImage3DDeformation &def)
{
    // Get the disp field from the def field
    NiftiImage3DTensor temp = def.deep_copy();
    reg_getDisplacementFromDeformation(temp.get_raw_nifti_sptr().get());
    temp.get_raw_nifti_sptr()->intent_p1 = DISP_FIELD;
    *this = temp.deep_copy();
}

void NiftiImage3DDisplacement::create_from_3D_image(const NiftiImage3D &image)
{
    this->NiftiImage3DTensor::create_from_3D_image(image);
    _nifti_image->intent_p1 = 1;
}

NiftiImage3DDeformation NiftiImage3DDisplacement::get_as_deformation_field(const NiftiImage3D &ref) const
{
    NiftiImage3DDeformation def;
    def.create_from_disp(*this);
    check_ref_and_def(ref,def);
    return def;
}
