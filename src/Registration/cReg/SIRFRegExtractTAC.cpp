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

#include "SIRFRegExtractTAC.h"
#include "SIRFRegMisc.h"

using namespace std;

void SIRFRegExtractTAC::update()
{
    // Check that all the necessary info has been set
    check_parameters();

    // Open the segmented image
    std::shared_ptr<nifti_image> segmented_sptr;
    SIRFRegMisc::open_nifti_image(segmented_sptr,_segmentation_filename);

    // Loop over each input image
    for (int i=0; i<_filenames.size(); i++) {

        // Open the file
        std::shared_ptr<nifti_image> im_sptr;
        SIRFRegMisc::open_nifti_image(im_sptr,_filenames[i]);


    }

}

void SIRFRegExtractTAC::check_parameters()
{
    // If anything is missing
    if (_filenames.size() == 0) {
        throw std::runtime_error("Input files have not been set.");}
    if (_segmentation_filename == "") {
        throw std::runtime_error("Segmentation has not been set."); }
    if (_VOIs.size() == 0) {
        throw std::runtime_error("VOIs have not been set."); }
}
