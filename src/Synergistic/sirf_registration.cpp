/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018-2019 University College London

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
\brief Perform registration with any types of SIRF image and any SIRF registration algorithm.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftyAladinSym.h"
#include "sirf/Reg/NiftyF3dSym.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/STIR/stir_data_containers.h"


using namespace sirf;

static std::shared_ptr<const ImageData> image_as_sptr(const std::string &filename, const std::string &engine = "Nifti")
{
    if      (strcmp(engine.c_str(), "Nifti") == 0)
        return std::make_shared<const NiftiImageData3D<float> >(filename);
    else if (strcmp(engine.c_str(), "STIR") == 0)
        return std::make_shared<const STIRImageData>(filename);
    else if (strcmp(engine.c_str(), "Gadgetron") == 0) {
        std::shared_ptr<GadgetronImageData> sptr_img(new GadgetronImagesVector);
		sptr_img->read(filename);
        return sptr_img;
    }
    else
        throw std::runtime_error("Synergistic_aladin: unknown image engine - " + engine + ".\n");
}

static std::shared_ptr<Registration<float> > algo_as_sptr(const std::string &algorithm)
{
    std::cout << "\nUsing " << algorithm << " registration algorithm...\n";
    if      (strcmp(algorithm.c_str(), "aladin") == 0)
        return std::make_shared<NiftyAladinSym<float> >();
    else if (strcmp(algorithm.c_str(), "f3d") == 0)
        return std::make_shared<NiftyF3dSym<float> >();
    else
        throw std::runtime_error("Synergistic_registration: unknown algorithm - " + algorithm + ".\n");
}

/// print usage
void print_usage()
{
    std::cout << "\n*** sirf_registration usage ***\n";

    // Required flags
    std::cout << "\n  Required flags:\n";
    std::cout << "    -algo:\t\tRegistration algorithm (aladin/f3d...)\n";
    std::cout << "    -ref:\t\treference image\n";
    std::cout << "    -flo:\t\tfloating image\n";
    std::cout << "    -par:\t\tparameter file\n";

    // Optional flags
    std::cout << "\n  Optional flags:\n";
    std::cout << "    -eng_ref:\t\tengine to open reference image (and reference mask if present) [Nifti|STIR|Gadgetron]\n";
    std::cout << "    -eng_flo:\t\tengine to open floating image (and floating mask if present) [Nifti|STIR|Gadgetron]\n";
    std::cout << "    -rmask:\t\tmask of reference image\n";
    std::cout << "    -fmask:\t\tmask of floating image\n";
    std::cout << "    -warped:\t\twarped image filename\n";
    std::cout << "    -TM_forward:\tforward transformation matrix\n";
    std::cout << "    -TM_inverse:\tinverse transformation matrix\n";
    std::cout << "    -disp_fwd:\tforward displacement field image\n";
    std::cout << "    -def_fwd:\tforward deformation field image\n";
    std::cout << "    -disp_inv:\tinverse displacement field image\n";
    std::cout << "    -def_inv:\tinverse deformation field image\n";
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


        // ------------------------------------------------ //
        // Required flags
        // ------------------------------------------------ //

        const std::string algo = argv[find_flag(unused_flags,argv,"-algo",true)+1];
        std::shared_ptr<Registration<float> > reg = algo_as_sptr(algo);

        // Get images
        int flag_ref = find_flag(unused_flags,argv,"-ref",true);
        int flag_eng_ref = find_flag(unused_flags,argv,"-eng_ref");
        std::shared_ptr<const ImageData> ref;
        if (flag_eng_ref==-1) {
            std::cout << "\nNo engine supplied for reference image, assuming Nifti.\n";
            ref = image_as_sptr(argv[flag_ref+1]);
        }
        else
            ref = image_as_sptr(argv[flag_ref+1],argv[flag_eng_ref+1]);

        int flag_flo = find_flag(unused_flags,argv,"-flo",true);
        int flag_eng_flo = find_flag(unused_flags,argv,"-eng_flo");
        std::shared_ptr<const ImageData> flo;
        if (flag_eng_flo==-1) {
            std::cout << "\nNo engine supplied for floating image, assuming Nifti.\n";
            flo = image_as_sptr(argv[flag_flo+1]);
        }
        else
            flo = image_as_sptr(argv[flag_flo+1],argv[flag_eng_flo+1]);

        // Set images
        reg->set_reference_image(ref);
        reg->set_floating_image(flo);

        int flag_par = find_flag(unused_flags,argv,"-par",true);
        reg->set_parameter_file(argv[flag_par+1]);


        // ------------------------------------------------ //
        // Optional flags
        // ------------------------------------------------ //

        // Masks
        int r_mask           = find_flag(unused_flags,argv,"-rmask");
        if (r_mask != -1)
            reg->set_reference_mask(image_as_sptr(argv[r_mask+1],argv[flag_eng_ref+1]));
        int f_mask           = find_flag(unused_flags,argv,"-fmask");
        if (f_mask != -1)
            reg->set_floating_mask(image_as_sptr(argv[f_mask+1],argv[flag_eng_flo+1]));

        // Warped image
        int flag_warped      = find_flag(unused_flags,argv,"-warped");

        // TMs
        int flag_TM_forward  = find_flag(unused_flags,argv,"-TM_forward");
        int flag_TM_inverse  = find_flag(unused_flags,argv,"-TM_inverse");

        // Forward/inverse disp/def field images
        int flag_disp_fwd = find_flag(unused_flags,argv,"-disp_fwd");
        int flag_disp_inv = find_flag(unused_flags,argv,"-disp_inv");
        int flag_def_fwd  = find_flag(unused_flags,argv,"-def_fwd");
        int flag_def_inv  = find_flag(unused_flags,argv,"-def_inv");


        // ------------------------------------------------ //
        // Unknown flags
        // ------------------------------------------------ //

        if (unused_flags.size() > 0) {
            std::cout << "\n\nThe following unknown flags were supplied:\n";
            for (unsigned i=0; i<unused_flags.size(); i++)
                std::cout << "\t" << argv[unused_flags[i]] << "\n";
            print_usage();
            return EXIT_FAILURE;
        }


        // ------------------------------------------------ //
        // Process
        // ------------------------------------------------ //

        reg->process();


        // ------------------------------------------------ //
        // Save results
        // ------------------------------------------------ //

        // Warped image
        if (flag_warped != -1)
            reg->get_output_sptr()->write(argv[flag_warped+1]);

        // TMs - only if rigid/affine
        if (strcmp(algo.c_str(),"aladin")==0) {
            if (flag_TM_forward != -1)
                std::dynamic_pointer_cast<NiftyAladinSym<float> >(reg)->get_transformation_matrix_forward_sptr()->write(argv[flag_TM_forward+1]);
            if (flag_TM_inverse != -1)
                std::dynamic_pointer_cast<NiftyAladinSym<float> >(reg)->get_transformation_matrix_inverse_sptr()->write(argv[flag_TM_inverse+1]);
        }

        // Forward disp field images
        if (flag_disp_fwd != -1)
            reg->get_displacement_field_forward_sptr()->write(argv[flag_disp_fwd+1]);

        // Forward def field images
        if (flag_def_fwd != -1)
            reg->get_deformation_field_forward_sptr()->write(argv[flag_def_fwd+1]);

        // Inverse disp field images
        if (flag_disp_inv != -1)
            reg->get_displacement_field_inverse_sptr()->write(argv[flag_disp_inv+1]);

        // Inverse def field images
        if (flag_def_inv != -1)
            reg->get_deformation_field_inverse_sptr()->write(argv[flag_def_inv+1]);
    }

    // If there was an error
    catch(const std::exception &error) {
        std::cerr << "\nError encountered:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
