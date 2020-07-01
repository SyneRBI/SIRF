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

#pragma once

#include <memory>

namespace sirf {

// Forward declarations
template<class dataType> class NiftiImageData;
template<class dataType> class NiftiImageData3DDeformation;

/*!
\ingroup Registration
\brief Class for converting control point grids to deformation field transformations.

\author Richard Brown
\author SyneRBI
*/
template<class dataType>
class ControlPointGridToDeformationConverter
{
public:

    /// CPG to DVF
    static std::shared_ptr<NiftiImageData3DDeformation<dataType> >
    forward(const NiftiImageData3DDeformation<dataType> &cpg_sptr,
            const NiftiImageData<dataType> &ref_sptr);

    /// CPG to DVF (in place)
    static void forward(NiftiImageData3DDeformation<dataType> &dvf,
                        const NiftiImageData3DDeformation<dataType> &cpg,
                        const NiftiImageData<dataType> &ref);

    /// DVF to CPG
    static std::shared_ptr<NiftiImageData3DDeformation<dataType> >
    backward(const NiftiImageData3DDeformation<dataType> &dvf,
             const NiftiImageData<dataType> &ref,
             float *spacingMillimeter);

    /// DVF to CPG (in place)
    static void backward(NiftiImageData3DDeformation<dataType> &cpg,
                         const NiftiImageData3DDeformation<dataType> &dvf,
                         const NiftiImageData<dataType> &ref,
                         float *spacingMillimeter);
};
}
