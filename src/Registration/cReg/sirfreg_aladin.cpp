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
\brief Perform niftyreg aladin registration

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFRegNiftyAladinSym.h"

using namespace sirf;

/// print usage
void print_usage()
{
    std::cout << "\n*** sirfreg_aladin usage ***\n";

    // Required flags
    std::cout << "\n  Required flags:\n";
    std::cout << "    -ref:\treference image\n";
    std::cout << "    -flo:\tfloating image\n";
    std::cout << "    -par:\tparameter file\n";

    // Optional flags
    std::cout << "\n  Optional flags:\n";
    std::cout << "    -warped:\twarped image filename\n";
    std::cout << "    -TM_forward:\tforward transformation matrix\n";
    std::cout << "    -TM_inverse:\tinverse transformation matrix\n";
    std::cout << "    -disp_4D:\t4D forward displacement field image\n";
    std::cout << "    -disp_3D:\t3D forward displacement field image\n";
    std::cout << "    -def_4D:\t4D forward deformation field image\n";
    std::cout << "    -def_3D:\t3D forward deformation field image\n";
}

/// find flag
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

        SIRFRegNiftyAladinSym<float> aladin;


        // ------------------------------------------------ //
        // Required flags
        // ------------------------------------------------ //

        int flag_ref = find_flag(unused_flags,argv,"-ref",true);
        NiftiImageData3D reference(argv[flag_ref+1]);
        aladin.set_reference_image(reference);

        int flag_flo = find_flag(unused_flags,argv,"-flo",true);
        NiftiImageData3D floating(argv[flag_flo+1]);
        aladin.set_floating_image(floating);

        int flag_par = find_flag(unused_flags,argv,"-par",true);
        aladin.set_parameter_file(argv[flag_par+1]);


        // ------------------------------------------------ //
        // Optional flags
        // ------------------------------------------------ //

        // Warped image
        int flag_warped   = find_flag(unused_flags,argv,"-warped");

        // TMs
        int flag_TM_forward  = find_flag(unused_flags,argv,"-TM_forward");
        int flag_TM_inverse  = find_flag(unused_flags,argv,"-TM_inverse");

        // Disp field images
        int flag_disp_4D = find_flag(unused_flags,argv,"-disp_4D");
        int flag_disp_3D = find_flag(unused_flags,argv,"-disp_3D");

        // Def field images
        int flag_def_4D  = find_flag(unused_flags,argv,"-def_4D");
        int flag_def_3D  = find_flag(unused_flags,argv,"-def_3D");


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
        // Process
        // ------------------------------------------------ //

        aladin.process();


        // ------------------------------------------------ //
        // Save results
        // ------------------------------------------------ //

        // Warped image
        if (flag_warped != -1)
            aladin.get_output().save_to_file(argv[flag_warped+1]);

        // TMs
        if (flag_TM_forward != -1)
            aladin.get_transformation_matrix_forward().save_to_file(argv[flag_TM_forward+1]);
        if (flag_TM_inverse != -1)
            aladin.get_transformation_matrix_inverse().save_to_file(argv[flag_TM_inverse+1]);

        // Disp field images
        if (flag_disp_4D != -1)
            aladin.get_displacement_field_forward().save_to_file(argv[flag_disp_4D+1]);
        if (flag_disp_3D != -1)
            aladin.get_displacement_field_forward().save_to_file_split_xyz_components(argv[flag_disp_3D+1]);

        // Def field images
        if (flag_def_4D != -1)
            aladin.get_deformation_field_forward().save_to_file(argv[flag_def_4D+1]);
        if (flag_def_3D != -1)
            aladin.get_deformation_field_forward().save_to_file_split_xyz_components(argv[flag_def_3D+1]);
    }

    // If there was an error
    catch(const std::exception &error) {
        std::cerr << "\nError encountered:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
