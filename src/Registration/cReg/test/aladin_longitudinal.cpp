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

#include "SIRFRegNiftyAladin.h"


#include "SIRFRegMisc.h"

using namespace std;

int aladin_longitudinal(string output_path)
{
    string SIRF_PATH     = getenv("SIRF_PATH");
    string examples_path = SIRF_PATH + "/data/examples/Registration";

    string reference_image_filename                = examples_path + "/test.nii.gz";
    string floating_image_filename                 = examples_path + "/test2.nii.gz";
    string parameter_file_aladin                   = examples_path + "/paramFiles/aladin.par";
    string warped_image_filename                   = output_path   + "/aladin_longitudinal_cplusplus";
    string transformation_matrix_filename          = output_path   + "/transformation_matrix_aladin_longitudinal_cplusplus.txt";
    string transformation_matrix_inversefilename   = output_path   + "/transformation_matrix_inverse_aladin_longitudinal_cplusplus.txt";

    // Run the test
    SIRFRegNiftyAladin<float> NA;
    NA.set_reference_image_filename      (        reference_image_filename       );
    NA.set_floating_image_filename       (        floating_image_filename        );
    NA.set_parameter_file                (         parameter_file_aladin         );
    NA.update();
    NA.save_warped_image                 (         warped_image_filename         );
    NA.save_transformation_matrix        (     transformation_matrix_filename    );
    NA.save_inverse_transformation_matrix( transformation_matrix_inversefilename );
}
