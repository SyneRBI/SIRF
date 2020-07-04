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
template<class dataType> class NiftiImageData3DBSpline;

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

    /// Constructor
    ControlPointGridToDeformationConverter();

    /// Set CPG spacing
    void set_cpg_spacing(const float spacing[3]);

    /// Set reference image for generating dvfs
    void set_reference_image(const NiftiImageData<dataType> &ref);

    /// CPG to DVF
    NiftiImageData3DDeformation<dataType> forward(const NiftiImageData3DBSpline<dataType> &cpg) const;

    /// DVF to CPG
    NiftiImageData3DBSpline<dataType> backward(const NiftiImageData3DDeformation<dataType> &dvf) const;

private:

    /// Check is set up
    void check_is_set_up() const;

    float _spacing[3];
    std::shared_ptr<NiftiImageData<dataType> > _template_ref_sptr;
};
}
