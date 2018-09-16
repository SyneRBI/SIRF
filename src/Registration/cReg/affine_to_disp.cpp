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
#include "SIRFRegMisc.h"
#include "SIRFImageDataDeformation.h"
#if NIFTYREG_VER_1_5
#include "_reg_globalTrans.h"
#endif

using namespace std;
using namespace sirf;

/// Print usage
void print_usage()
{
    cout << "\n*** affine_to_disp usage ***\n";

    // Required flags
    cout << "\n  Required flags:\n";
    cout << "    -TM:\taffine transformation matrix\n";
    cout << "    -ref:\treference image\n";

    // Optional flags
    cout << "\n  Optional flags (but at least one is required):\n";
    cout << "    -disp_4D:\t4D forward displacement field image\n";
    cout << "    -disp_3D:\t3D forward displacement field image\n";
    cout << "    -def_4D:\t4D forward deformation field image\n";
    cout << "    -def_3D:\t3D forward deformation field image\n";
}

/// Find flag
int find_flag(vector<int> &unused_flags, char* argv[], string arg, bool required=false)
{
    for (int i=0; i<unused_flags.size(); i++) {
        if (!strcmp(argv[unused_flags[i]], arg.c_str())) {
            int flag = unused_flags[i];
            unused_flags.erase(unused_flags.begin() + i);
            return flag;
        }
    }

    if (!required)
        return -1;

    cout << "\nRequired flag \"" << arg << "\" not found. Exiting...\n";
    print_usage();
    exit(EXIT_FAILURE);
}

/// main
int main(int argc, char* argv[])
{

    try {
        // Create a list of all unused flags
        vector<int> unused_flags;
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
        string TM_filename = argv[flag_TM+1];

        int flag_ref = find_flag(unused_flags,argv,"-ref",true);
        string ref_filename = argv[flag_ref+1];


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
            cout << "\nRequired at least one output type. Exiting...\n";
            print_usage();
            return EXIT_FAILURE;
        }


        // ------------------------------------------------ //
        // Unknown flags
        // ------------------------------------------------ //

        if (unused_flags.size() > 0) {
            cout << "\n\nThe following unknown flags were supplied:\n";
            for (int i=0; i<unused_flags.size(); i++)
                cout << "\t" << argv[unused_flags[i]] << "\n";
        }
        cout << "\n";


        // ------------------------------------------------ //
        // Do the conversion
        // ------------------------------------------------ //

        // Open the transformation matrix
        mat44 TM;
        SIRFRegMisc::open_transformation_matrix(TM,TM_filename);

        // Create images
        std::shared_ptr<nifti_image> ref_sptr, cpp_sptr;
        SIRFImageDataDeformation def, disp;

        // Open reference image
        SIRFRegMisc::open_nifti_image(ref_sptr, ref_filename);

        // Get the deformation field image
#if NIFTYREG_VER_1_5
        def.create_from_3D_image(ref_sptr);
        reg_affine_getDeformationField(&TM, def.get_image_as_nifti().get());
#elif NIFTYREG_VER_1_3
        SIRFRegMisc::get_cpp_from_transformation_matrix(cpp_sptr, TM, ref_sptr);
        SIRFRegMisc::get_def_from_cpp(def_sptr,cpp_sptr, ref_sptr);
#endif

        // Get the displacement fields from the def
        disp = def;
        SIRFRegMisc::convert_from_def_to_disp(disp);

        // If they want to save the deformation field images
        if (flag_def_4D != -1)
            def.save_to_file(argv[flag_def_4D+1],false,   "4D deformation field");
        if (flag_def_3D != -1)
            def.save_to_file(argv[flag_def_3D+1],true,    "3D deformation field");
        // If they want to save the displacement field images
        if (flag_disp_4D != -1)
            disp.save_to_file(argv[flag_disp_4D+1],false, "4D displacement field");
        if (flag_disp_3D != -1)
            disp.save_to_file(argv[flag_disp_3D+1],true,  "3D displacement field");


    // If there was an error
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
