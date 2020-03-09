/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 University College London.

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

#include "sirf/common/ImageData.h"

using namespace sirf;

void ImageData::reorient(const VoxelisedGeometricalInfo3D &)
{
    throw std::runtime_error("ImageData::reorient not yet implemented for your image type.");
}

bool ImageData::can_reorient(const VoxelisedGeometricalInfo3D &geom_1, const VoxelisedGeometricalInfo3D &geom_2, const bool throw_error)
{
    // If size and spacing match, return true
    if (geom_1.get_size() == geom_2.get_size() && geom_1.get_spacing() == geom_2.get_spacing())
        return true;
    // Else (and error desired), print error
    if (throw_error)
        throw std::runtime_error("ImageData::can_reorient: num voxels or spacing do not match.");
    // Else, return false
    return false;
}