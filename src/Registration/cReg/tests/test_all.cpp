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

#include <iostream>
#include "activity_correction.cpp"
#include "aladin_longitudinal.cpp"
#include "aladin_multimodal.cpp"
#include "f3d_mouse.cpp"
#include "resample.cpp"
#include "weighted_mean.cpp"

using namespace std;

int main(int argc, char* argv[])
{

    try {

    	string output_path = argv[0];
    	output_path = output_path.substr(0, output_path.find_last_of('/'));
        
        activity_correction(output_path);
		aladin_longitudinal(output_path);
		aladin_multimodal(output_path);
		f3d_mouse(output_path);
		resample(output_path);
		weighted_mean(output_path);

    // If there was an error
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
