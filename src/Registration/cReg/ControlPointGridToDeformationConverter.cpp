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

using namespace sirf;

template<class dataType>
std::shared_ptr<NiftiImageData3DDeformation<dataType> >
ControlPointGridToDeformationConverter<dataType>::
forward(const NiftiImageData3DDeformation<dataType> &cpg,
        const NiftiImageData<dataType> &ref)
{
    std::shared_ptr<NiftiImageData3DDeformation<dataType> > dvf_sptr =
            std::make_shared<NiftiImageData3DDeformation<dataType> >();
    forward(*dvf_sptr, cpg, ref);
    return dvf_sptr;
}

template<class dataType>
void
ControlPointGridToDeformationConverter<dataType>::
forward(NiftiImageData3DDeformation<dataType> &dvf,
        const NiftiImageData3DDeformation<dataType> &cpg,
        const NiftiImageData<dataType> &ref
        )
{
    dvf.create_from_cpp(cpg, ref);
}

template<class dataType>
std::shared_ptr<NiftiImageData3DDeformation<dataType> >
ControlPointGridToDeformationConverter<dataType>::
backward(const NiftiImageData3DDeformation<dataType> &dvf,
         const NiftiImageData<dataType> &ref,
         float *spacingMillimeter)
{
    std::shared_ptr<NiftiImageData3DDeformation<dataType> > cpg_sptr =
            std::make_shared<NiftiImageData3DDeformation<dataType> >();
    dvf.create_cpp(*cpg_sptr, ref, spacingMillimeter);
    return cpg_sptr;
}

template<class dataType>
void
ControlPointGridToDeformationConverter<dataType>::
backward(NiftiImageData3DDeformation<dataType> &cpg,
         const NiftiImageData3DDeformation<dataType> &dvf,
         const NiftiImageData<dataType> &ref,
         float *spacingMillimeter)
{
    dvf.create_cpp(cpg, ref, spacingMillimeter);
}

namespace sirf {
template class ControlPointGridToDeformationConverter<float>;
}
