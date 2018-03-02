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
\brief Test for Registration_Nifty_f3d (non-linear registration with NiftyReg)

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFRegNiftyF3d.h"

using namespace std;

int f3d_mouse(string output_path)
{
    string SIRF_PATH     = getenv("SIRF_PATH");
    string examples_path = SIRF_PATH + "/data/examples/Registration";
    
    string reference_image_filename = examples_path + "/mouseFixed.nii.gz";
    string floating_image_filename  = examples_path + "/mouseMoving.nii.gz";
    string parameter_file_f3d       = examples_path + "/paramFiles/f3d.par";
    string warped_image_filename    = output_path   + "/f3d_mice_cplusplus";

    // Run the test
    SIRFRegNiftyF3d<float> NF;

    NF.set_reference_image_filename(reference_image_filename);
    NF.set_floating_image_filename (floating_image_filename );
    NF.set_parameter_file          (   parameter_file_f3d   );
    NF.set_reference_time_point    (            1           );
    NF.set_floating_time_point     (            1           );
    NF.update();
    NF.save_warped_image           ( warped_image_filename  );
}
