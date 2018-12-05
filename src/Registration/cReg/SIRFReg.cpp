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
\brief Base class for all SIRF registration.

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFReg.h"
#include <nifti1_io.h>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace sirf;

template<class dataType>
void SIRFReg<dataType>::check_parameters() const
{
    // If anything is missing
    if (_parameter_filename.empty())
        throw std::runtime_error("Parameter file has not been set.");
    if (!_floating_image_sptr->is_initialised())
        throw std::runtime_error("Floating image has not been set.");
    if (!_reference_image_sptr->is_initialised())
        throw std::runtime_error("Reference image has not been set.");
}

template<class dataType>
void SIRFReg<dataType>::set_parameter(const std::string &par, const std::string &arg1, const std::string &arg2)
{
    _extra_params.push_back(par);
    _extra_params.push_back(arg1);
    _extra_params.push_back(arg2);
}

template<class dataType>
const std::shared_ptr<const NiftiImageData3D<dataType> > SIRFReg<dataType>::get_output() const
{
    return std::const_pointer_cast<const NiftiImageData3D<dataType> >(this->_warped_image_sptr);
}

template<class dataType>
const std::shared_ptr<const NiftiImageData3DDeformation<dataType> > SIRFReg<dataType>::get_deformation_field_forward() const
{
    NiftiImageData3DDeformation<dataType> def_fwd;
    def_fwd.create_from_disp(*this->_disp_image_forward_sptr);
    return std::make_shared<const NiftiImageData3DDeformation<dataType> >(def_fwd);
}

template<class dataType>
const std::shared_ptr<const NiftiImageData3DDeformation<dataType> > SIRFReg<dataType>::get_deformation_field_inverse() const
{
    NiftiImageData3DDeformation<dataType> def_inv;
    def_inv.create_from_disp(*this->_disp_image_inverse_sptr);
    return std::make_shared<const NiftiImageData3DDeformation<dataType> >(def_inv);
}

template<class dataType>
const std::shared_ptr<const NiftiImageData3DDisplacement<dataType> > SIRFReg<dataType>::get_displacement_field_forward() const
{
    return std::const_pointer_cast<const NiftiImageData3DDisplacement<dataType> >(this->_disp_image_forward_sptr);
}
template<class dataType>
const std::shared_ptr<const NiftiImageData3DDisplacement<dataType> > SIRFReg<dataType>::get_displacement_field_inverse() const
{
    return std::const_pointer_cast<const NiftiImageData3DDisplacement<dataType> >(this->_disp_image_inverse_sptr);
}

namespace sirf {
template class SIRFReg<float>;
}
