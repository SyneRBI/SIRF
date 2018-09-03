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

#include "SIRFRegTransformation.h"
#include "SIRFRegMisc.h"
#include <_reg_globalTrans.h>

using namespace std;

SIRFImageDataDeformation SIRFRegTransformationAffine::get_as_deformation_field(const SIRFImageData &ref)
{
    SIRFImageDataDeformation def;
    def.create_from_3D_image(ref);
    reg_affine_getDeformationField(&_tm, def.get_image_as_nifti().get());
    def.get_image_as_nifti()->intent_p1 = DEF_FIELD;
    return def;
}

SIRFImageDataDeformation SIRFRegTransformationDisplacement::get_as_deformation_field(const SIRFImageData &)
{
    SIRFImageDataDeformation def;
    def = _disp.deep_copy();
    SIRFRegMisc::convert_from_disp_to_def(def);
    return def;
}

SIRFImageDataDeformation SIRFRegTransformationDeformation::get_as_deformation_field(const SIRFImageData &)
{
    return _def;
}

SIRFRegTransformationAffine SIRFRegTransformationAffine::deep_copy() const
{
    SIRFRegTransformationAffine temp(_tm);
    return temp;
}

SIRFRegTransformationDisplacement SIRFRegTransformationDisplacement::deep_copy() const
{
    SIRFRegTransformationDisplacement temp(_disp);
    return temp;
}

SIRFRegTransformationDeformation SIRFRegTransformationDeformation::deep_copy() const
{
    SIRFRegTransformationDeformation temp(_def);
    return temp;
}
