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
\brief Class for converting control point grids to deformation field transformations.

\author Richard Brown
\author SyneRBI
*/

#include "sirf/Reg/ControlPointGridToDeformationConverter.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/NiftiImageData3DBSpline.h"
#include "sirf/NiftyMoMo/BSplineTransformation.h"

using namespace sirf;

template<class dataType>
ControlPointGridToDeformationConverter<dataType>::
ControlPointGridToDeformationConverter()
{
    for (unsigned i=0; i<3; ++i)
        _spacing[i] = std::numeric_limits<float>::quiet_NaN();
}

template<class dataType>
void
ControlPointGridToDeformationConverter<dataType>::
set_cpg_spacing(const float spacing[3])
{
    for (unsigned i=0; i<3; ++i)
        _spacing[i] = spacing[i];
}

template<class dataType>
void
ControlPointGridToDeformationConverter<dataType>::
set_reference_image(const NiftiImageData<dataType> &ref)
{
    _template_ref_sptr = ref.clone();
}

template<class dataType>
NiftiImageData3DDeformation<dataType>
ControlPointGridToDeformationConverter<dataType>::
forward(const NiftiImageData3DBSpline<dataType> &cpg) const
{
    check_is_set_up();
//    NiftiImageData3DDeformation<float> dvf;
//    dvf.create_from_cpp(cpg, *_template_ref_sptr);
//    return dvf;
    return cpg.get_as_deformation_field(*_template_ref_sptr);
}

template<class dataType>
NiftiImageData3DBSpline<dataType>
ControlPointGridToDeformationConverter<dataType>::
backward(const NiftiImageData3DDeformation<dataType> &dvf) const
{
    check_is_set_up();
    // not marked const, so copy
    float spacing_nonconst[3] = {_spacing[0], _spacing[1], _spacing[2]};
    // Get raw nifti_image from reference image
    nifti_image *ref_ptr = _template_ref_sptr->get_raw_nifti_sptr().get();
    // Create the NiftyMoMo bspline transformation class
    NiftyMoMo::BSplineTransformation bspline(ref_ptr, 1, spacing_nonconst);
    // Get cpg_ptr
    nifti_image *cpg_ptr = bspline.GetTransformationAsImage();
    // Convert DVF to CPG
    std::shared_ptr<NiftiImageData3DDeformation<dataType> > dvf_sptr = dvf.clone();
    cpg_ptr->data = bspline.GetDVFGradientWRTTransformationParameters(dvf_sptr->get_raw_nifti_sptr().get());
    cpg_ptr->intent_p1 = SPLINE_VEL_GRID;
    return NiftiImageData3DBSpline<dataType>(*cpg_ptr);
}

template<class dataType>
void ControlPointGridToDeformationConverter<dataType>::
check_is_set_up() const
{
    // Has spacing been set?
    for (unsigned i=0; i<3; ++i)
        if (std::isnan(_spacing[i]))
            throw std::runtime_error("ControlPointGridToDeformationConverter: Set CPG spacing.");

    // Has template deformation been set?
    if (!_template_ref_sptr)
        throw std::runtime_error("ControlPointGridToDeformationConverter: Set template DVF.");
}

namespace sirf {
template class ControlPointGridToDeformationConverter<float>;
}
