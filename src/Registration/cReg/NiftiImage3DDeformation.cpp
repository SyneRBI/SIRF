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
#include "NiftiImage3DDeformation.h"
#include "NiftiImage3DDisplacement.h"
#include "SIRFRegMisc.h"
#include <_reg_globalTrans.h>
#include <sstream>
#include <_reg_localTrans.h>

using namespace std;
using namespace sirf;

void NiftiImage3DDeformation::create_from_disp(const NiftiImage3DDisplacement &disp)
{
    // Get the def field from the disp field
    NiftiImage3DTensor temp = disp.deep_copy();
    reg_getDeformationFromDisplacement(temp.get_raw_nifti_sptr().get());
    temp.get_raw_nifti_sptr()->intent_p1 = DEF_FIELD;
    *this = temp.deep_copy();
}

void NiftiImage3DDeformation::create_from_3D_image(const NiftiImage3D &image)
{
    this->NiftiImage3DTensor::create_from_3D_image(image);
    //_nifti_image->intent_p1 = 0; not necessary. 0 by default
}

void NiftiImage3DDeformation::create_from_cpp(NiftiImage3DTensor &cpp, const NiftiImage3D &ref)
{
    this->create_from_3D_image(ref);

    reg_spline_getDeformationField(cpp.get_raw_nifti_sptr().get(),
                                   _nifti_image.get(),
                                   NULL,
                                   false, //composition
                                   true // bspline
                                   );
}


NiftiImage3DDeformation NiftiImage3DDeformation::get_as_deformation_field(const NiftiImage3D &ref) const
{
    check_ref_and_def(ref,*this);
    return this->deep_copy();
}

NiftiImage3DDeformation NiftiImage3DDeformation::compose_single_deformation(const vector<SIRFRegTransformation*> &transformations, const NiftiImage3D &ref)
{
    if (transformations.size() == 0)
        throw runtime_error("NiftiImage3DDeformation::compose_single_deformation no transformations given.");

    NiftiImage3DDeformation def = transformations.at(0)->get_as_deformation_field(ref).deep_copy();

    for (unsigned i=1; i<transformations.size(); ++i) {
        NiftiImage3DDeformation temp = transformations.at(i)->get_as_deformation_field(ref);
        reg_defField_compose(temp.get_raw_nifti_sptr().get(),def.get_raw_nifti_sptr().get(),nullptr);
    }
    return def;
}

NiftiImage3DDeformation NiftiImage3DDeformation::compose_single_deformation(const vector<shared_ptr<SIRFRegTransformation> > &transformations, const NiftiImage3D &ref)
{
    vector<SIRFRegTransformation*> vec;
    for (unsigned i=0; i<transformations.size(); ++i)
        vec.push_back(transformations.at(i).get());
    return compose_single_deformation(vec, ref);
}
