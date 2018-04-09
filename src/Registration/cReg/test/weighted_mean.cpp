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

#include "SIRFRegImageWeightedMean.h"
#include <iostream>

using namespace std;

int weighted_mean(string output_path)
{
    string SIRF_PATH     = getenv("SIRF_PATH");
    string examples_path = SIRF_PATH + "/data/examples/Registration";
    
    string im1    = examples_path + "/weighted_mean/regis_recon_gate1.nii";
    string im2    = examples_path + "/weighted_mean/regis_recon_gate2.nii";
    string im3    = examples_path + "/weighted_mean/regis_recon_gate3.nii";
    string im4    = examples_path + "/weighted_mean/regis_recon_gate4.nii";
    string output = output_path   + "/weighted_mean";

    // Run the test
    SIRFRegImageWeightedMean weighted_mean;
    weighted_mean.add_image( im1, 0.2 );
    weighted_mean.add_image( im2, 0.2 );
    weighted_mean.add_image( im3, 0.2 );
    weighted_mean.add_image( im4, 0.2 );
    weighted_mean.update();
    weighted_mean.save_image_to_file(output);
}
