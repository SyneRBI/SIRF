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
\brief Split or join tensor images (i.e., deformation or displacement fields)

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include "NiftiImage3DTensor.h"
#include "NiftiImage3D.h"

using namespace std;
using namespace sirf;

void print_usage()
{
    cout << "\nUsage: sirfreg_tensor_split_join --join/split filename_4D filename_x filename_y filename z\n";
}

enum JoinOrSplit{join,split};

/// main
int main(int argc, char* argv[])
{
    try {
        if (argc != 5) {
            print_usage();
            EXIT_SUCCESS;
        }

        // Are we splitting or joining images?
        JoinOrSplit mode = join;
        if (strcmp(argv[1], "--split") == 0)
            mode = split;
        else if (strcmp(argv[1], "--join") == 0)
            mode = join;
        else {
            print_usage();
            EXIT_SUCCESS;
        }

        // Get filenames
        string filename_4D = argv[2];
        string filename_x  = argv[3];
        string filename_y  = argv[4];
        string filename_z  = argv[5];

        // If we're joining images
        if (mode == join) {
            std::cout << "\nDoing join.\n";
            NiftiImage3D x(filename_x);
            NiftiImage3D y(filename_x);
            NiftiImage3D z(filename_x);
            NiftiImage3DTensor tensor(x,y,z);
            tensor.save_to_file(filename_4D);
        }
        // If we're splitting
        else {
            std::cout << "\nDoing split.\n";
            NiftiImage3DTensor tensor(filename_4D);
            tensor.save_to_file_split_xyz_components(filename_x,filename_y,filename_z);
        }

    // If there was an error
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
