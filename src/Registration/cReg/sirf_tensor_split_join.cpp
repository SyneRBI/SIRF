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
\brief Split or join tensor images (i.e., deformation or displacement fields)

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include "sirf/Reg/NiftiImageData3DTensor.h"
#include "sirf/Reg/NiftiImageData3D.h"

using namespace sirf;

void print_usage()
{
    std::cout << "\nUsage: sirf_tensor_split_join --join/split filename_4D filename_x filename_y filename_z\n";
}

enum JoinOrSplit{join,split};

/// main
int main(int argc, char* argv[])
{
    try {
        if (argc != 6) {
            print_usage();
            return EXIT_FAILURE;
        }

        // Are we splitting or joining images?
        JoinOrSplit mode = join;
        if (strcmp(argv[1], "--split") == 0)
            mode = split;
        else if (strcmp(argv[1], "--join") == 0)
            mode = join;
        else {
            print_usage();
            return EXIT_FAILURE;
        }

        // Get filenames
        std::string filename_4D = argv[2];
        std::string filename_x  = argv[3];
        std::string filename_y  = argv[4];
        std::string filename_z  = argv[5];

        // If we're joining images
        if (mode == join) {
            std::cout << "\nDoing join.\n";
            NiftiImageData3D<float> x(filename_x);
            NiftiImageData3D<float> y(filename_x);
            NiftiImageData3D<float> z(filename_x);
            NiftiImageData3DTensor<float> tensor(x,y,z);
            tensor.write(filename_4D);
        }
        // If we're splitting
        else {
            std::cout << "\nDoing split.\n";
            NiftiImageData3DTensor<float> tensor(filename_4D);
            tensor.write_split_xyz_components(filename_x,filename_y,filename_z);
        }

    // If there was an error
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
