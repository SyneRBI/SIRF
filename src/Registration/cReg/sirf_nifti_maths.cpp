/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2019 University College London

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
\brief Crop a nifti image

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include <vector>
#include <sirf/Reg/NiftiImageData.h>

using namespace sirf;

/// Print usage
void print_usage()
{
    std::cout << "\n\n\n*** Usage: sirf_nifti_maths output_filename --<type> [type_param] input_filename_1 [input_filename_2][-h/--help]***\n\n";
    std::cout << "\t--<type> describes the maths type and can be one of:\n";
    std::cout << "\t\t --add_scalar, --mul_scalar, --sub_scalar, --div_scalar. A value must follow immediately after these flags and only 1 input image should be given.\n";
    std::cout << "\t\t --add, --sub. A second image must be given for these options.\n";

    std::cout << "\nExample usage:\n";
    std::cout << "1. sirf_nifti_maths out --add_scalar 5 input.nii\n";
    std::cout << "2. sirf_nifti_maths out --add input1.nii input2.nii\n";
}

/// main
int main(int argc, char* argv[])
{

    try {

        // Check for help
        for (int i=1; i<argc; ++i) {
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                print_usage();
                return EXIT_SUCCESS;
            }
        }

        // Always 5 inputs
        if (argc != 5) {
            print_usage();
            return EXIT_FAILURE;
        }

        // Output image
        const std::string output_filename = argv[1];

        const std::string operation = argv[2];
        bool is_scalar;
        if      (operation.compare("--add_scalar") == 0) is_scalar = true;
        else if (operation.compare("--sub_scalar") == 0) is_scalar = true;
        else if (operation.compare("--div_scalar") == 0) is_scalar = true;
        else if (operation.compare("--mul_scalar") == 0) is_scalar = true;
        else if (operation.compare("--add")        == 0) is_scalar = false;
        else if (operation.compare("--sub")        == 0) is_scalar = false;
        else {
            std::cout << "\n\nsirf_nifti_maths: unknown operation\n\n";
            print_usage();
            return EXIT_FAILURE;
        }

        if (is_scalar) {
            float val = std::stof(argv[3]);
            const std::string input_filename = argv[4];
            NiftiImageData<float> im(input_filename);
            if      (operation.compare("--add_scalar") == 0)
                im += val;
            else if (operation.compare("--sub_scalar") == 0)
                im -= val;
            else if (operation.compare("--div_scalar") == 0)
                im /= val;
            else if (operation.compare("--mul_scalar") == 0)
                im *= val;
            im.write(output_filename);
        }
        else {
            const std::string input_filename_1 = argv[3];
            const std::string input_filename_2 = argv[4];
            NiftiImageData<float> im1(input_filename_1);
            NiftiImageData<float> im2(input_filename_2);
            if      (operation.compare("--add") == 0)
                im1 += im2;
            else if (operation.compare("--sub") == 0)
                im1 -= im2;
            im1.write(output_filename);
        }

    // If there was an error
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
