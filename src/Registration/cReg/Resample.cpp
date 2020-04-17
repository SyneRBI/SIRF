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
\brief Abstract resampling base class

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/Resample.h"

using namespace sirf;

/// Set reference image
template<class dataType>
void Resample<dataType>::set_reference_image(const std::shared_ptr<const ImageData> reference_image_sptr)
{
    _reference_image_sptr = reference_image_sptr;
    _need_to_set_up = true;
    _need_to_set_up_forward = true;
    _need_to_set_up_adjoint = true;
}

/// Set floating image
template<class dataType>
void Resample<dataType>::set_floating_image(const std::shared_ptr<const ImageData> floating_image_sptr)
{
    _floating_image_sptr = floating_image_sptr;
    _need_to_set_up = true;
    _need_to_set_up_forward = true;
    _need_to_set_up_adjoint = true;
}

template<class dataType>
void Resample<dataType>::add_transformation(const std::shared_ptr<const Transformation<dataType> > transformation_sptr)
{
    _transformations.push_back(transformation_sptr);
    this->_need_to_set_up = true;
    _need_to_set_up_forward = true;
    _need_to_set_up_adjoint = true;
}

template<class dataType>
void Resample<dataType>::set_interpolation_type(const enum InterpolationType type)
{
    _interpolation_type = type;
    this->_need_to_set_up = true;
    _need_to_set_up_forward = true;
    _need_to_set_up_adjoint = true;
}

template<class dataType>
void Resample<dataType>::check_parameters()
{
    // If anything is missing
    if (!_reference_image_sptr)
        throw std::runtime_error("Reference image has not been set.");
    if (!_floating_image_sptr)
        throw std::runtime_error("Floating image has not been set.");
    if (_interpolation_type == NOTSET)
        throw std::runtime_error("Interpolation type has not been set.");
}

template<class dataType>
std::shared_ptr<ImageData> Resample<dataType>::backward(const std::shared_ptr<const ImageData> input_sptr)
{
    return adjoint(input_sptr);
}

template<class dataType>
void Resample<dataType>::backward(std::shared_ptr<ImageData> output_sptr, const std::shared_ptr<const ImageData> input_sptr)
{
    adjoint(output_sptr, input_sptr);
}

namespace sirf {
template class Resample<float>;
}

