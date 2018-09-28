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
#include <sstream>

using namespace std;
using namespace sirf;

void SIRFRegTransformation::check_ref_and_def(const NiftiImage3D &ref, const NiftiImage3DDeformation &def) const
{
    // Check image size of ref matches def
    int ref_dims[8], def_dims[8];
    ref.get_dimensions(ref_dims);
    def.get_dimensions(def_dims);

    bool all_ok = true;
    for (int i=1; i<=3; ++i)
        if (ref_dims[i] != def_dims[i])
            all_ok = false;

    if (!all_ok) {
        stringstream ss;
        ss << "Deformation field image should contain same number of x, y and z voxels.\n";
        ss << "Reference: ";
        for (int i=1; i<=3; ++i) ss << ref_dims[i] << " ";
        ss << "\nDeformation: ";
        for (int i=1; i<=3; ++i) ss << def_dims[i] << " ";
        throw std::runtime_error(ss.str());
    }
}
