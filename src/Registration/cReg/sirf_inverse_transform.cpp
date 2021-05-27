/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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
\author SyneRBI
*/

#include <sirf/Reg/NiftiImageData.h>
#include <sirf/Reg/AffineTransformation.h>
#include <sirf/Reg/NiftiImageData3DDisplacement.h>

using namespace sirf;

/// Print usage
static void print_usage_and_exit(const char * const program_name, const int exit_status)
{
    std::cout << "\n\n\n*** Usage: " << program_name << " out_filename in_filename <type> [--flo <filename>][--use_vtk][-h/--help]***\n\n";
    std::cout << "\t<type>             can be \"aff\", \"disp\" or \"def\".\n";
    std::cout << "\t[--flo <filename>] if output type is non-rigid, optionally supply a floating image on which transformation will be resampled.\n";
    std::cout << "\t[--use_vtk]        if SIRF was built with VTK, --useVTK can be used when inverting a displacement field image.\n";

    std::cout << "\nExample usage:\n";
    std::cout << "1. sirf_inverse_transform out input.txt aff\n";
    std::cout << "2. sirf_inverse_transform out input.nii disp --flo floating.nii\n";
    std::cout << "3. sirf_inverse_transform out input.nii def --useVTK\n";
    exit(exit_status);
}

/// main
int main(int argc, char* argv[])
{

    try {
        const char * const program_name = argv[0];

        // Check for help request
        for (int i=1; i<argc; ++i)
            if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0)
                print_usage_and_exit(program_name, EXIT_SUCCESS);

        // Check for all compulsory arguments
        if (argc<5)
            print_usage_and_exit(program_name, EXIT_FAILURE);
        // Output image
        const std::string output_filename = argv[1];
        const std::string input_filename  = argv[2];
        const std::string type            = argv[3];

        // skip past compulsory arguments
        argc-=4;
        argv+=4;

        // Set default value for optional arguments
        std::shared_ptr<const NiftiImageData<float> > floating_sptr;
        bool use_vtk = false;

        // Loop over remaining input
        while (argc>0 && argv[0][0]=='-') {
            if (strcmp(argv[0], "--flo")==0) {
                std::string floating_filename = argv[1];
                floating_sptr = std::make_shared<const NiftiImageData<float> >(floating_filename);
                argc-=2; argv+=2;
            }
            else if (strcmp(argv[0], "--use_vtk")==0) {
                use_vtk = true;
                argc-=1; argv+=1;
            }
            else {
                std::cerr << "Unknown option '" << argv[0] <<"'\n";
                print_usage_and_exit(program_name, EXIT_FAILURE);
            }
        }

        // Affine
        if (type.compare("aff")==0) {
            const AffineTransformation<float> in_trans(input_filename);
            in_trans.get_inverse().write(output_filename);
        }

        // Deformation
        else if (type.compare("def")==0) {
            const NiftiImageData3DDeformation<float> in_trans(input_filename);
            in_trans.get_inverse(floating_sptr,use_vtk)->write(output_filename);
        }

        // Displacement
        else if (type.compare("disp")==0) {
            const NiftiImageData3DDisplacement<float> in_trans(input_filename);
            in_trans.get_inverse(floating_sptr,use_vtk)->write(output_filename);
        }

        // Unknown transformation type
        else
            throw std::runtime_error("Unknown input type.");

    // If there was an error
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
