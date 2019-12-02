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
\brief Flip or mirror image

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include <sirf/Reg/NiftiImageData.h>

using namespace sirf;

/// Print usage
void print_usage(const char * program_name)
{
    std::cout << "\n\n\n*** Usage: " << program_name << " [-h] output_filename input_filename [--flip <axis_to_flip> | --mirror <axis_to_mirror>]***\n\n";
    std::cout << "\t(By flip we mean rotation about 180 degrees, whereas by mirror we mean to switch the handedness.)\n";
}

/// main
int main(int argc, char* argv[])
{

    try {

        for (int i=1; i<argc; ++i) {
            if (strcmp(argv[i],"-h")==0) {
                print_usage(argv[0]);
                return EXIT_SUCCESS;
            }
        }

        if (argc != 5) {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        const std::string output_filename = argv[1];
        const std::string input_filename = argv[2];
        const unsigned axis = unsigned(std::stoul(argv[4]));

        // Read image
        NiftiImageData<float> im(input_filename);

        // If mirroring the image
        if (strcmp(argv[3],"--flip")==0) {
            im.flip_along_axis(axis);
        }
        else if (strcmp(argv[3],"--mirror")==0) {
            im.mirror_along_axis(axis);
        }
        else {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        // Save to file
        im.write(output_filename);


    // If there was an error
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
