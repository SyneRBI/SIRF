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
\brief Convert affine transformation matrix to displacement or deformation field image(s)

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include <vector>
#include <NiftiImage.h>

using namespace sirf;

/// Print usage
void print_usage()
{
    std::cout << "\n\n\n*** Usage: sirfreg_crop_image output_filename input_filename [at least one optional parameter] ***\n\n";
    std::cout << "Where optional parameters are:\n";
    std::cout << "\t[--x_min <val>]\n";
    std::cout << "\t[--x_max <val>]\n";
    std::cout << "\t[--y_min <val>]\n";
    std::cout << "\t[--y_max <val>]\n";
    std::cout << "\t[--z_min <val>]\n";
    std::cout << "\t[--z_max <val>]\n";
    std::cout << "\t[--t_min <val>]\n";
    std::cout << "\t[--t_max <val>]\n";
    std::cout << "\t[--u_min <val>]\n";
    std::cout << "\t[--u_max <val>]\n";
    std::cout << "\t[--v_min <val>]\n";
    std::cout << "\t[--v_max <val>]\n";
    std::cout << "\t[--w_min <val>]\n";
    std::cout << "\t[--w_max <val>]\n";
    std::cout << "When optional parameters are not supplied, the image will remain unchanged for that aspect.\n\n";
}

/// main
int main(int argc, char* argv[])
{

    try {

        // We should always more than 3 inputs, and always odd number
        if (argc < 4 && argc % 2 == 1) {
            print_usage();
            return EXIT_SUCCESS;
        }

        const std::string input_filename  = argv[2];
        const std::string output_filename = argv[1];

        // Ignore first 3 arguments
        argv+=3;
        argc-=3;

        // Open the image
        NiftiImage im(input_filename);

        // Get the dims of the nifti image
        const int *dims = im.get_dimensions();

        // Set the default values; if the user doesn't specify otherwise, the
        // image will be uncropped in that direction
        int min_index[7], max_index[7];
        for (int i=0; i<7; ++i) {
            min_index[i] = 0;
            max_index[i] = dims[i+1] - 1;
        }

        // Read new min and max indices
        while (argc>0 && argv[0][0]=='-') {

            // x
            if      (strcmp(argv[0], "--x_min")==0)
                min_index[0] = atoi(argv[1]);
            else if (strcmp(argv[0], "--x_max")==0)
                max_index[0] = atoi(argv[1]);

            // y
            else if (strcmp(argv[0], "--y_min")==0)
                min_index[1] = atoi(argv[1]);
            else if (strcmp(argv[0], "--y_max")==0)
                max_index[1] = atoi(argv[1]);

            // z
            else if (strcmp(argv[0], "--z_min")==0)
                min_index[2] = atoi(argv[1]);
            else if (strcmp(argv[0], "--z_max")==0)
                max_index[2] = atoi(argv[1]);

            // t
            else if (strcmp(argv[0], "--t_min")==0)
                min_index[3] = atoi(argv[1]);
            else if (strcmp(argv[0], "--t_max")==0)
                max_index[3] = atoi(argv[1]);

            // u
            else if (strcmp(argv[0], "--u_min")==0)
                min_index[4] = atoi(argv[1]);
            else if (strcmp(argv[0], "--u_max")==0)
                max_index[4] = atoi(argv[1]);

            // v
            else if (strcmp(argv[0], "--v_min")==0)
                min_index[5] = atoi(argv[1]);
            else if (strcmp(argv[0], "--v_max")==0)
                max_index[5] = atoi(argv[1]);

            // w
            else if (strcmp(argv[0], "--w_min")==0)
                min_index[6] = atoi(argv[1]);
            else if (strcmp(argv[0], "--w_max")==0)
                max_index[6] = atoi(argv[1]);

            // Unknown
            else {
                std::cerr << "Unknown option '" << argv[0] <<"'\n";
                return EXIT_FAILURE;
            }

            argc-=2;
            argv+=2;
        }

        std::cout << "\nDimensions of input image = ( ";
        for (int i=1; i<8; ++i)
            std::cout << dims[i] << " ";
        std::cout << ")\nDesired minimum index     = ( ";
        for (int i=0; i<7; ++i)
            std::cout << min_index[i] << " ";
        std::cout << ")\nDesired maximum index     = ( ";
        for (int i=0; i<7; ++i)
            std::cout << max_index[i] << " ";
        std::cout << ")\n\n";

        // Crop
        im.crop(min_index,max_index);

        // Save output
        im.save_to_file(output_filename);

    // If there was an error
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
