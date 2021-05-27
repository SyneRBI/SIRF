/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2019 University College London

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
\brief Base class for transformations.

\author Richard Brown
\author SyneRBI
*/

#include "sirf/Reg/Transformation.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include <sstream>

using namespace sirf;

template<class dataType>
void Transformation<dataType>::check_ref_and_def(const NiftiImageData<dataType> &ref, const NiftiImageData3DDeformation<dataType> &def)
{
    // Check image size of ref matches def
    const int *ref_dims = ref.get_dimensions();
    const int *def_dims = def.get_dimensions();

    bool all_ok = true;
    for (int i=1; i<=3; ++i)
        if (ref_dims[i] != def_dims[i])
            all_ok = false;

    if (!all_ok) {
        std::stringstream ss;
        ss << "Deformation field image should contain same number of x, y and z voxels.\n";
        ss << "Reference: ";
        for (int i=1; i<=3; ++i) ss << ref_dims[i] << " ";
        ss << "\nDeformation: ";
        for (int i=1; i<=3; ++i) ss << def_dims[i] << " ";
        throw std::runtime_error(ss.str());
    }
}

namespace sirf {
template class Transformation<float>;
}
