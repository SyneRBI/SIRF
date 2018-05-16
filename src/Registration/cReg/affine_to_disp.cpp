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
#include "SIRFRegMisc.h"

using namespace std;

int main(int argc, char* argv[])
{

    try {
        if (argc<4 || argc>5) {
            std::cerr << "\n\tUsage: " << argv[0] << " <1> <2> <3> <4> <5>\n";
            std::cerr << "\t<1>: Transformation matrix filename\n";
            std::cerr << "\t<2>: reference image filename\n";
            std::cerr << "\t<3>: Output image filename\n";
            std::cerr << "\t<4>: Split xyz components of output image (0/1. optional, defaults to 1)\n";
            std::cerr << "\t<5>: Flip for STIR (0/1. optional, defaults to 1)\n";
            return EXIT_FAILURE;
        }

        string TM_filename        = argv[1];
        string ref_filename       = argv[2];
        string output_filename    = argv[3];
        bool   split_xyz          = true;
        bool   flip_for_STIR      = true;
        if (argc>4) split_xyz     = atoi(argv[4]);
        if (argc>5) flip_for_STIR = atoi(argv[5]);

        // Open the transformation matrix
        std::shared_ptr<mat44> TM_sptr;
        SIRFRegMisc::open_transformation_matrix(TM_sptr,TM_filename);

        // Create images
        shared_ptr<nifti_image> ref_sptr, cpp_sptr, def_sptr, disp_sptr;

        // Open reference image
        SIRFRegMisc::open_nifti_image(ref_sptr, ref_filename);

        // Create control point image
        SIRFRegMisc::get_cpp_from_transformation_matrix(cpp_sptr, TM_sptr, ref_sptr);

        // Get deformation fields from cpp
        SIRFRegMisc::get_def_from_cpp(def_sptr,cpp_sptr, ref_sptr);

        // Get the displacement fields from the def
        SIRFRegMisc::get_disp_from_def(disp_sptr,def_sptr);

        // Write the displacement field image
        SIRFRegMisc::save_multicomponent_nifti_image(disp_sptr,output_filename,split_xyz,flip_for_STIR);

    // If there was an error
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
