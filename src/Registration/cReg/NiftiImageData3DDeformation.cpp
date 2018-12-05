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

#include "NiftiImageData3DDeformation.h"
#include "SIRFRegTransformation.h"
#include "NiftiImageData3DDisplacement.h"
#include <_reg_globalTrans.h>
#include <sstream>
#include <_reg_localTrans.h>

using namespace sirf;

template<class dataType>
void NiftiImageData3DDeformation<dataType>::create_from_disp(const NiftiImageData3DDisplacement<dataType> &disp)
{
    // Get the def field from the disp field
    NiftiImageData3DTensor<dataType> temp = disp.deep_copy();
    reg_getDeformationFromDisplacement(temp.get_raw_nifti_sptr().get());
    temp.get_raw_nifti_sptr()->intent_p1 = DEF_FIELD;
    *this = temp.deep_copy();
}

template<class dataType>
void NiftiImageData3DDeformation<dataType>::create_from_3D_image(const NiftiImageData3D<dataType> &image)
{
    this->NiftiImageData3DTensor<dataType>::create_from_3D_image(image);
    //_nifti_image->intent_p1 = 0; not necessary. 0 by default
}

template<class dataType>
void NiftiImageData3DDeformation<dataType>::create_from_cpp(NiftiImageData3DTensor<dataType> &cpp, const NiftiImageData3D<dataType> &ref)
{
    this->create_from_3D_image(ref);

    reg_spline_getDeformationField(cpp.get_raw_nifti_sptr().get(),
                                   this->_nifti_image.get(),
                                   NULL,
                                   false, //composition
                                   true // bspline
                                   );
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DDeformation<dataType>::get_as_deformation_field(const NiftiImageData3D<dataType> &ref) const
{
    this->check_ref_and_def(ref,*this);
    return this->deep_copy();
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DDeformation<dataType>::compose_single_deformation(const std::vector<const SIRFRegTransformation<dataType>*> &transformations, const NiftiImageData3D<dataType> &ref)
{
    if (transformations.size() == 0)
        throw std::runtime_error("NiftiImageData3DDeformation::compose_single_deformation no transformations given.");

    NiftiImageData3DDeformation def = transformations.at(0)->get_as_deformation_field(ref).deep_copy();

    for (unsigned i=1; i<transformations.size(); ++i) {
        NiftiImageData3DDeformation temp = transformations.at(i)->get_as_deformation_field(ref);
        reg_defField_compose(temp.get_raw_nifti_sptr().get(),def.get_raw_nifti_sptr().get(),nullptr);
    }
    return def;
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DDeformation<dataType>::compose_single_deformation(const std::vector<std::shared_ptr<const SIRFRegTransformation<dataType> > > &transformations, const NiftiImageData3D<dataType> &ref)
{
    std::vector<const SIRFRegTransformation<dataType>*> vec;
    for (unsigned i=0; i<transformations.size(); ++i)
        vec.push_back(transformations.at(i).get());
    return compose_single_deformation(vec, ref);
}

template<class dataType>
std::shared_ptr<SIRFRegTransformation<dataType> > NiftiImageData3DDeformation<dataType>::get_clone_sptr() const
{
    return std::shared_ptr<NiftiImageData3DDeformation<dataType> >
            (new NiftiImageData3DDeformation<dataType>(this->deep_copy()));
}

namespace sirf {
template class NiftiImageData3DDeformation<float>;
}
