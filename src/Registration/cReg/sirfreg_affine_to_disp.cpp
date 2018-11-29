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
#include "NiftiImageData3D.h"
#include "NiftiImageData3DDisplacement.h"
#include "NiftiImageData3DDeformation.h"
#include "SIRFRegAffineTransformation.h"

using namespace sirf;

/// Print usage
void print_usage()
{
    std::cout << "\n*** sirfreg_affine_to_disp usage ***\n";

    // Required flags
    std::cout << "\n  Required flags:\n";
    std::cout << "    -TM:\taffine transformation matrix\n";
    std::cout << "    -ref:\treference image\n";

    // Optional flags
    std::cout << "\n  Optional flags (but at least one is required):\n";
    std::cout << "    -disp_4D:\t4D forward displacement field image\n";
    std::cout << "    -disp_3D:\t3D forward displacement field image\n";
    std::cout << "    -def_4D:\t4D forward deformation field image\n";
    std::cout << "    -def_3D:\t3D forward deformation field image\n";

    std::cout << "If 3D format is required, use boost format for the filename (e.g., output_%s.nii).\n\n\n";
}

/// Find flag
int find_flag(std::vector<int> &unused_flags, char* argv[], std::string arg, bool required=false)
{
    for (unsigned i=0; i<unused_flags.size(); i++) {
        if (!strcmp(argv[unused_flags[i]], arg.c_str())) {
            int flag = unused_flags[i];
            unused_flags.erase(unused_flags.begin() + i);
            return flag;
        }
    }

    if (!required)
        return -1;

    std::cout << "\nRequired flag \"" << arg << "\" not found. Exiting...\n";
    print_usage();
    exit(EXIT_FAILURE);
}

/// main
int main(int argc, char* argv[])
{

    try {
        // Create a list of all unused flags
        std::vector<int> unused_flags;
        for (int i=1; i<argc; i+=2)
            unused_flags.push_back(i);

        // Print help if desired
        if (find_flag(unused_flags,argv,"-h") + find_flag(unused_flags,argv,"-help") > -2) {
            print_usage();
            return EXIT_SUCCESS;
        }



        // ------------------------------------------------ //
        // Required flags
        // ------------------------------------------------ //

        int flag_TM = find_flag(unused_flags,argv,"-TM",true);
        std::string TM_filename = argv[flag_TM+1];

        int flag_ref = find_flag(unused_flags,argv,"-ref",true);
        std::string ref_filename = argv[flag_ref+1];


        // ------------------------------------------------ //
        // Optional flags
        // ------------------------------------------------ //

        // Disp field images
        int flag_disp_4D = find_flag(unused_flags,argv,"-disp_4D");
        int flag_disp_3D = find_flag(unused_flags,argv,"-disp_3D");

        // Def field images
        int flag_def_4D  = find_flag(unused_flags,argv,"-def_4D");
        int flag_def_3D  = find_flag(unused_flags,argv,"-def_3D");

        // If all are blank nothing to do
        if (flag_def_3D == -1 &&
                flag_def_4D  == -1 &&
                flag_disp_3D == -1 &&
                flag_disp_4D == -1) {
            std::cout << "\nRequired at least one output type. Exiting...\n";
            print_usage();
            return EXIT_FAILURE;
        }


        // ------------------------------------------------ //
        // Unknown flags
        // ------------------------------------------------ //

        if (unused_flags.size() > 0) {
            std::cout << "\n\nThe following unknown flags were supplied:\n";
            for (unsigned i=0; i<unused_flags.size(); i++)
                std::cout << "\t" << argv[unused_flags[i]] << "\n";
        }
        std::cout << "\n";


        // ------------------------------------------------ //
        // Do the conversion
        // ------------------------------------------------ //

        // Open the transformation matrix
        SIRFRegAffineTransformation TM(TM_filename);

        // Create images

        // Open reference image
        NiftiImageData3D ref(ref_filename);

        // Get the deformation field image
        NiftiImageData3DDeformation def = TM.get_as_deformation_field(ref);

        // Get the displacement fields from the def
        NiftiImageData3DDisplacement disp;
        disp.create_from_def(def);

        // If they want to save the deformation field images
        if (flag_def_4D != -1)
            def.save_to_file(argv[flag_def_4D+1]);
        if (flag_def_3D != -1)
            def.save_to_file_split_xyz_components(argv[flag_def_3D+1]);
        // If they want to save the displacement field images
        if (flag_disp_4D != -1)
            disp.save_to_file(argv[flag_disp_4D+1]);
        if (flag_disp_3D != -1)
            disp.save_to_file_split_xyz_components(argv[flag_disp_3D+1]);


    // If there was an error
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
