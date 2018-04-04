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

#include "SIRFRegNiftyResample.h"
#include <iostream>


#include "SIRFRegMisc.h"

using namespace std;

int resample(string output_path)
{
    cout << "\n========================================================\n";
    cout << "    TESTING RESAMPLING";
    cout << "\n========================================================\n";

    string SIRF_PATH     = getenv("SIRF_PATH");
    string examples_path = SIRF_PATH + "/data/examples/Registration";

    string reference = examples_path + "/test.nii.gz";
    string floating  = examples_path + "/test2.nii.gz";
    string matrix    = examples_path + "/transformation_matrix.txt";
    string output    = output_path   + "/resampled_image_CPLUSPLUS.nii";

    // Run the test
    SIRFRegNiftyResample resample;
    resample.set_reference_image_filename       ( reference );
    resample.set_floating_image_filename        ( floating  );
    resample.add_transformation_matrix_filename (  matrix   );
    resample.set_interpolation_type_to_cubic_spline();
    resample.update();

    resample.save_resampled_image               (  output   );

    cout << "\n========================================================\n";
    cout << "    SUCCESSFULLY COMPLETED TESTING RESAMPLING";
    cout << "\n========================================================\n";

    return 0;
}
