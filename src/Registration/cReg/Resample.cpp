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
\brief bstract resampling base class

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/cReg/Resample.h"

using namespace sirf;

template<class dataType>
void SIRFRegResample<dataType>::add_transformation(const std::shared_ptr<const SIRFRegTransformation<dataType> > transformation_sptr)
{
    _transformations.push_back(transformation_sptr);
}

template<class dataType>
void SIRFRegResample<dataType>::check_parameters()
{
    // If anything is missing
    if (!_reference_image_sptr)
        throw std::runtime_error("Reference image has not been set.");
    if (!_floating_image_sptr)
        throw std::runtime_error("Floating image has not been set.");
    if (_interpolation_type == NOTSET)
        throw std::runtime_error("Interpolation type has not been set.");
}

namespace sirf {
template class SIRFRegResample<float>;
}

