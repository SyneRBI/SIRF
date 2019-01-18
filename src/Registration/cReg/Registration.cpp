/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2019 University College London

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
\brief Base class for all SIRF registration.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/cReg/Registration.h"
#include <stdexcept>

using namespace sirf;

template<class dataType>
void Registration<dataType>::check_parameters() const
{
    // If anything is missing
    if (!_floating_image_sptr)
        throw std::runtime_error("Floating image has not been set.");
    if (!_reference_image_sptr)
        throw std::runtime_error("Reference image has not been set.");
}

template<class dataType>
void Registration<dataType>::set_parameter(const std::string &par, const std::string &arg1, const std::string &arg2)
{
    _extra_params.push_back(par);
    _extra_params.push_back(arg1);
    _extra_params.push_back(arg2);
}

namespace sirf {
template class Registration<float>;
}
