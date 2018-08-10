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
#include "SIRFRegMisc.h"
#include <nifti1_io.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <boost/filesystem.hpp>

using namespace std;

void SIRFReg::check_parameters()
{
    // If anything is missing
    if (_parameter_filename == "") {
        throw std::runtime_error("Parameter file has not been set.");}
    if (!_floating_image.is_initialised()) {
        throw std::runtime_error("Floating image has not been set."); }
    if (!_reference_image.is_initialised()) {
        throw std::runtime_error("Reference image has not been set."); }
}

void SIRFReg::save_warped_image(const string filename) const
{
    _warped_image.save_to_file(filename);
}

void SIRFReg::save_deformation_field_fwrd_image(const std::string &filename, const bool &split_xyz)
{
    _def_image_fwrd.save_to_file(filename,split_xyz,"fowrard deformation");
}

void SIRFReg::save_deformation_field_back_image(const std::string &filename, const bool &split_xyz)
{
    _def_image_back.save_to_file(filename,split_xyz,"backwards deformation");
}

void SIRFReg::save_displacement_field_fwrd_image(const std::string &filename, const bool &split_xyz)
{
    _disp_image_fwrd.save_to_file(filename,split_xyz,"forward displacement");
}

void SIRFReg::save_displacement_field_back_image(const std::string &filename, const bool &split_xyz)
{
    _disp_image_back.save_to_file(filename,split_xyz,"backwards displacement");
}
