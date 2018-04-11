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
#include <_reg_localTransformation.h>
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
    if (!_floating_image_sptr && _floating_image_filename == "") {
        throw std::runtime_error("Floating image has not been set."); }
    if (!_reference_image_sptr && _reference_image_filename == "") {
        throw std::runtime_error("Reference image has not been set."); }
}

void SIRFReg::save_warped_image(const string filename) const
{
    SIRFRegMisc::save_nifti_image(_warped_image_sptr,filename);
}

void SIRFReg::save_displacement_field_fwrd_image(const std::string &filename, const bool &split_xyz, const bool &flip_for_stir)
{
    save_displacement_field_image(_disp_image_fwrd_sptr,filename,split_xyz,flip_for_stir);
}

void SIRFReg::save_displacement_field_back_image(const std::string &filename, const bool &split_xyz, const bool &flip_for_stir)
{
    save_displacement_field_image(_disp_image_back_sptr,filename,split_xyz,flip_for_stir);
}

void SIRFReg::save_displacement_field_image(const std::shared_ptr<nifti_image> &im_sptr, const std::string &filename, const bool &split_xyz, const bool &flip_for_stir)
{
    // Check that the disp image exists
    if (!im_sptr)
        throw std::runtime_error("Displacement field image not available. Have you run the registration?");

    // Check that filename isn't blank
    if (filename == "")
        throw std::runtime_error("Error, cannot write displacement field image to file because filename is blank.");

    cout << "\nSaving displacement field image to file (" << filename << ")..." << flush;

    shared_ptr<nifti_image> im_to_save_sptr;

    // If the user wants to save it as it is
    if (!flip_for_stir)
        im_to_save_sptr = im_sptr;

    // But if they want to output for stir, need to flip the z-axis
    else {
        im_to_save_sptr = make_shared<nifti_image>();
        // Deep copy the displacement image
        SIRFRegMisc::copy_nifti_image(im_to_save_sptr,im_sptr);
        // Flip the z-axis
        SIRFRegMisc::flip_multicomponent_image(im_to_save_sptr,2);
    }

    // Whether anything has been flipped or not, save the result...

    // If the user wants it saved as multicomponent image
    if (!split_xyz)
        SIRFRegMisc::save_nifti_image(im_to_save_sptr,filename);

    // If the user wants the multicomponent image split into 3 separate images.
    else
        SIRFRegMisc::save_split_multicomponent_nifti_image(im_to_save_sptr,filename);

    cout << "Done.\n";
}
