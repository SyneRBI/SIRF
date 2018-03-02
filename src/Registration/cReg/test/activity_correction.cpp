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
\brief Test for creation of a weighted-mean image

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFRegActivityCorrect.h"
#include <iostream>

using namespace std;

int activity_correction(string output_path)
{
    string SIRF_PATH     = getenv("SIRF_PATH");
    string examples_path = SIRF_PATH + "/data/examples/Registration";

    string input  = examples_path + "/test.nii.gz";
    string output = output_path   + "/activity_correction_CPLUSPLUS.nii";

    // Run the test
    SIRFRegActivityCorrect act_corr;
    act_corr.set_initial_activity(267000000);
    act_corr.set_half_life(6586.2);
    act_corr.set_input_image_filename(input);
    act_corr.set_start(48);
    act_corr.set_stop(96);
    act_corr.update();
    act_corr.save_output(output);
}
