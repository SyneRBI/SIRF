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

#include <sirf/Reg/NiftiImageData.h>
#include <sirf/Reg/AffineTransformation.h>
#include <sirf/Reg/NiftiImageData3DDisplacement.h>

using namespace sirf;

/// Print usage
void print_usage()
{
    std::cout << "\n\n\n*** Usage: sirf_inverse_transform out_filename <type> in_filename <type> [floating_filename][-h/--help]***\n\n";
    std::cout << "\t<type> can be \"aff\", \"disp\" or \"def\".\n";
    std::cout << "\tA floating image is required when the output is a non-rigid transformation, and should describe the floating image where the inverted transformation is defined.\n";

    std::cout << "\nExample usage:\n";
    std::cout << "1. sirf_inverse_transform out aff input.txt aff\n";
    std::cout << "2. sirf_inverse_transform out disp input.nii def floating.nii\n";
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

        // arg should be 5 or 6
        if (argc < 5 || argc > 6) {
            print_usage();
            return EXIT_FAILURE;
        }

        // Output image
        const std::string output_filename   = argv[1];
        const std::string output_type       = argv[2];
        const std::string input_filename    = argv[3];
        const std::string input_type        = argv[4];

        // Check for existing types
        if (!(input_type.compare("def") == 0 || input_type.compare("disp") == 0 || input_type.compare("aff") == 0))
            throw std::runtime_error("Unknown input type.");
        if (!(output_type.compare("def") == 0 || output_type.compare("disp") == 0 || output_type.compare("aff") == 0))
            throw std::runtime_error("Unknown output type.");

        bool input_is_affine  = input_type.compare("aff") == 0;
        bool output_is_affine = input_type.compare("aff") == 0;

        // If output is non-rigid, need a floating image
        std::shared_ptr<const NiftiImageData<float> > floating_sptr;
        if (!input_is_affine) {
            if (argc != 6) {
                print_usage();
                return EXIT_FAILURE;
            }
            const std::string floating_filename = argv[5];
            floating_sptr = std::make_shared<const NiftiImageData<float> >(floating_filename);
        }

        std::shared_ptr<const NiftiImageData3DDeformation<float> > deformation_sptr;

        // If input is affine
        if (input_is_affine) {
            const std::shared_ptr<const AffineTransformation<float> > in_trans_sptr = std::make_shared<const AffineTransformation<float> >(input_filename);
            // If output is affine
            if (output_is_affine) {
                in_trans_sptr->get_inverse().write(output_filename);
                return EXIT_SUCCESS;
            }
            deformation_sptr = std::make_shared<const NiftiImageData3DDeformation<float> >(
                        in_trans_sptr->get_inverse().get_as_deformation_field(*floating_sptr));
        }
        // If non-rigid
        else {

            // If in transformation is non-rigid, out can't be affine.
            if (output_is_affine)
                throw std::runtime_error("Input transformation is non-rigid. Output can't be affine.");

            // If disp, read it and convert to deformation
            if (input_type.compare("disp") == 0) {
                const std::shared_ptr<const NiftiImageData3DDisplacement<float> > disp_sptr = std::make_shared<const NiftiImageData3DDisplacement<float> >(input_filename);
                deformation_sptr = std::make_shared<const NiftiImageData3DDeformation<float> >(*disp_sptr);
            }
            else
                deformation_sptr = std::make_shared<const NiftiImageData3DDeformation<float> >(input_filename);
            // Inverse
            deformation_sptr = deformation_sptr->get_inverse(floating_sptr);
        }

        // For the output, convert
        if (output_type.compare("def") == 0)
            deformation_sptr->write(output_filename);
        else
            std::make_shared<const NiftiImageData3DDisplacement<float> >(*deformation_sptr)->write(output_filename);

    // If there was an error
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
