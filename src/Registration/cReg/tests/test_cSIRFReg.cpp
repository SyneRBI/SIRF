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
\brief Test for creation of a weighted-mean image

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include "SIRFRegNiftyAladinSym.h"
#include "SIRFRegNiftyF3dSym.h"
#include "SIRFRegNiftyResample.h"
#include "SIRFImageData.h"
#include "SIRFRegImageWeightedMean.h"
#include "stir_data_containers.h"
#include <memory>

using namespace std;
using namespace sirf;

int main(int argc, char* argv[])
{

    try {

    // Paths
    string SIRF_PATH;
    if (argc==1)
        SIRF_PATH = getenv("SIRF_PATH");
    else
        SIRF_PATH = argv[1];
    string examples_path = SIRF_PATH + "/data/examples/Registration";
    string output_path   = "results/cplusplus_";

    // Input filenames
    string ref_aladin_filename      = examples_path + "/test.nii.gz";
    string flo_aladin_filename      = examples_path + "/test2.nii.gz";
    string ref_f3d_filename         = examples_path + "/mouseFixed.nii.gz";
    string flo_f3d_filename         = examples_path + "/mouseMoving.nii.gz";
    string parameter_file_aladin    = examples_path + "/paramFiles/aladin.par";
    string parameter_file_f3d       = examples_path + "/paramFiles/f3d.par";
    string matrix                   = examples_path + "/transformation_matrix.txt";
    string stir_nifti               = examples_path + "/nifti_created_by_stir.nii";

    // Output filenames
    string img_data                 = output_path + "save_SIRFImageData";
    string img_data_def_not_split   = output_path + "save_SIRFImageDataDeformation_not_split";
    string img_data_def_split       = output_path + "save_SIRFImageDataDeformation_split";
    string aladin_warped            = output_path   + "aladin_warped";
    string f3d_warped               = output_path   + "f3d_warped";
    string TM_fwrd                  = output_path   + "TM_fwrd.txt";
    string TM_back                  = output_path   + "TM_back.txt";
    string aladin_def_fwrd          = output_path   + "aladin_def_fwrd";
    string aladin_def_back          = output_path   + "aladin_def_back";
    string aladin_disp_fwrd         = output_path   + "aladin_disp_fwrd";
    string aladin_disp_back         = output_path   + "aladin_disp_back";
    string f3d_disp_fwrd            = output_path   + "f3d_disp_fwrd";
    string f3d_disp_back            = output_path   + "f3d_disp_back";
    string f3d_def_fwrd             = output_path   + "f3d_def_fwrd";
    string f3d_def_back             = output_path   + "f3d_def_back";
    string rigid_resample           = output_path   + "rigid_resample";
    string nonrigid_resample_disp   = output_path   + "nonrigid_resample_disp";
    string nonrigid_resample_def    = output_path   + "nonrigid_resample_def";
    string output_weighted_mean     = output_path   + "weighted_mean";
    string output_weighted_mean_def = output_path   + "weighted_mean_def";

    string output_stir_nifti        = output_path   + "stir_nifti.nii";

    SIRFImageData ref_aladin( ref_aladin_filename );
    SIRFImageData flo_aladin( flo_aladin_filename );
    SIRFImageData ref_f3d   (   ref_f3d_filename  );
    SIRFImageData flo_f3d   (   flo_f3d_filename  );
    SIRFImageData nifti     (      stir_nifti     );

    float required_percentage_accuracy = 1.F;

    bool try_misc_functions = true;
    bool try_sirfimagedata = true;
    bool try_sirfimagedatadeformation = true;
    bool try_niftyaladin = true;
    bool try_niftyf3d = true;
    bool try_transformations = true;
    bool try_resample = true;
    bool try_weighted_mean = true;
    bool try_stir_to_sirfreg = true;

    if (try_misc_functions) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting misc functions test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // do nifti images match?
        if (!SIRFRegMisc::do_nifti_images_match(ref_aladin, ref_aladin, required_percentage_accuracy))
            throw runtime_error("Images don't match, but they should.");
        cerr << "\nThe following images intentionally do not match.\n";
        if (SIRFRegMisc::do_nifti_images_match(ref_aladin, flo_aladin, required_percentage_accuracy))
            throw runtime_error("Images don't match, but they should.");

        // dump from filename
        SIRFRegMisc::dump_nifti_info(ref_aladin_filename);
        // dump from SIRFImageData
        SIRFRegMisc::dump_nifti_info(ref_aladin);
        // dump from multiple images
        std::vector<SIRFImageData> vec;
        vec.push_back(ref_aladin);
        vec.push_back(flo_aladin);
        vec.push_back(nifti);
        SIRFRegMisc::dump_nifti_info(vec);
        // dump from SIRFImageDataDeformation
        SIRFImageDataDeformation deform;
        deform.create_from_3D_image(ref_aladin);
        SIRFRegMisc::dump_nifti_info(deform);

        // identity matrix
        mat44 tm_iden = SIRFRegMisc::get_identity_matrix();

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished misc functions test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_sirfimagedata) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting SIRFImageData test...\n";
        cout << "//------------------------------------------------------------------------ //\n";


        // default constructor
        SIRFImageData a;

        // Read from file
        SIRFImageData b(ref_aladin_filename);

        // Construct from stir pSTIR.ImageData
        sirf::PETImageData stir(stir_nifti);
        SIRFImageData c(stir);

        // Save to file
        c.save_to_file(img_data);

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw runtime_error("SIRFImageData fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw runtime_error("SIRFImageData fill()/get_min() failed.");

        // Copy data to pSTIRImageData
        stir.fill(3.);
        c.copy_data_to(stir);
        if (fabs(stir.data_sptr()->find_max() - 100) > 1.e-5F)
            throw runtime_error("SIRFImageData copy_data_to stir ImageData failed.");

        // Deep copy
        SIRFImageData d = c.deep_copy();
        if (&d == &c)
            throw runtime_error("SIRFImageData deep_copy failed.");
        if (!SIRFRegMisc::do_nifti_images_match(d, c, required_percentage_accuracy))
            throw runtime_error("SIRFImageData deep_copy failed.");

        // Addition
        SIRFImageData e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw runtime_error("SIRFImageData __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw runtime_error("SIRFImageData __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw runtime_error("SIRFImageData get_sum() failed.");

        // Dimensions
        int f[8];
        int g[8] = {3, 285, 285, 127, 1, 1, 1, 1};
        e.get_dimensions(f);
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw runtime_error("SIRFImageData get_dimensions() failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished SIRFImageData test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_sirfimagedatadeformation) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting SIRFImageDataDeformation test...\n";
        cout << "//------------------------------------------------------------------------ //\n";


        // Create SIRFImageDataDeformation from SIRFImageData
        SIRFImageDataDeformation b;
        b.create_from_3D_image(ref_aladin);

        // Save to file
        b.save_to_file(img_data_def_not_split, false);
        b.save_to_file(img_data_def_split, true);

        // Constructor from file
        SIRFImageDataDeformation c(img_data_def_not_split + ".nii");

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw runtime_error("SIRFImageDataDeformation fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw runtime_error("SIRFImageDataDeformation fill()/get_min() failed.");

        // Deep copy
        SIRFImageDataDeformation d = c.deep_copy();
        if (&d == &c)
            throw runtime_error("SIRFImageDataDeformation deep_copy failed.");
        if (!SIRFRegMisc::do_nifti_images_match(d, c, required_percentage_accuracy))
            throw runtime_error("SIRFImageDataDeformation deep_copy failed.");

        // Addition
        SIRFImageData e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw runtime_error("SIRFImageDataDeformation __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw runtime_error("SIRFImageDataDeformation __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw runtime_error("SIRFImageDataDeformation get_sum() failed.");

        // Dimensions
        int f[8];
        int g[8] = {5, 64, 64, 64, 1, 3, 1, 1};
        e.get_dimensions(f);
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw runtime_error("SIRFImageDataDeformation get_dimensions() failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished SIRFImageDataDeformation test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    SIRFRegNiftyAladinSym<float> NA;
    if (try_niftyaladin) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting Nifty aladin test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        NA.set_reference_image               (           ref_aladin          );
        NA.set_floating_image                (           flo_aladin          );
        NA.set_parameter_file                (      parameter_file_aladin    );
        NA.update();
        NA.get_output().save_to_file         (         aladin_warped         );
        NA.save_transformation_matrix_fwrd   (             TM_fwrd           );
        NA.save_transformation_matrix_back   (             TM_back           );
        NA.get_displacement_field_fwrd().save_to_file( aladin_disp_fwrd, false       );
        NA.get_displacement_field_back().save_to_file( aladin_disp_back, true        );
        NA.get_deformation_field_fwrd().save_to_file( aladin_def_fwrd,  false       );
        NA.get_deformation_field_back().save_to_file( aladin_def_back,  true        );

        // Get outputs
        SIRFImageData warped = NA.get_output();
        SIRFImageDataDeformation def_fwrd  = NA.get_deformation_field_fwrd();
        SIRFImageDataDeformation def_back  = NA.get_deformation_field_back();
        SIRFImageDataDeformation disp_fwrd = NA.get_displacement_field_fwrd();
        SIRFImageDataDeformation disp_back = NA.get_displacement_field_back();

        // Fwrd TM
        mat44 fwrd_tm = NA.get_transformation_matrix_fwrd();
        SIRFRegMisc::print_mat44(fwrd_tm);

        // Back TM
        mat44 back_tm = NA.get_transformation_matrix_back();
        SIRFRegMisc::print_mat44(back_tm);

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished Nifty aladin test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }


    if (try_niftyf3d)
    {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting Nifty f3d test..\n";
        cout << "//------------------------------------------------------------------------ //\n";

        SIRFRegNiftyF3dSym<float> NF;
        NF.set_reference_image               (          ref_f3d           );
        NF.set_floating_image                (          flo_f3d           );
        NF.set_parameter_file                (     parameter_file_f3d     );
        NF.set_reference_time_point          (             1              );
        NF.set_floating_time_point           (             1              );
        NF.update();
        NF.get_output().save_to_file         (         f3d_warped         );
        NF.get_deformation_field_fwrd().save_to_file ( f3d_def_fwrd,  true        );
        NF.get_deformation_field_back().save_to_file ( f3d_def_back,  false       );
        NF.get_displacement_field_fwrd().save_to_file( f3d_disp_fwrd, true        );
        NF.get_displacement_field_back().save_to_file( f3d_disp_back, false       );

        // Get outputs
        SIRFImageData warped = NF.get_output();
        SIRFImageDataDeformation def_fwrd  = NF.get_deformation_field_fwrd();
        SIRFImageDataDeformation def_back  = NF.get_deformation_field_back();
        SIRFImageDataDeformation disp_fwrd = NF.get_displacement_field_fwrd();
        SIRFImageDataDeformation disp_back = NF.get_displacement_field_back();

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished Nifty f3d test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_transformations) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting transformations test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // Affine
        cout << "\nTesting affine...\n";
        SIRFRegTransformationAffine a1;
        SIRFRegTransformationAffine a2(TM_fwrd);
        SIRFRegTransformationAffine a3(NA.get_transformation_matrix_fwrd());

        // Displacement
        cout << "\nTesting displacement...\n";
        SIRFRegTransformationDisplacement b1;
        SIRFRegTransformationDisplacement b2(aladin_disp_fwrd + ".nii");
        SIRFRegTransformationDisplacement b3(NA.get_displacement_field_fwrd());

        // Deformation
        cout << "\nTesting deformation...\n";
        SIRFRegTransformationDeformation c1;
        SIRFRegTransformationDeformation c2(aladin_def_fwrd + ".nii");
        SIRFRegTransformationDeformation c3(NA.get_deformation_field_fwrd());

        // Get as deformations
        SIRFImageDataDeformation a_def = a3.get_as_deformation_field(ref_aladin);
        SIRFImageDataDeformation b_def = b3.get_as_deformation_field(ref_aladin);
        SIRFImageDataDeformation c_def = c3.get_as_deformation_field(ref_aladin);
        if (!SIRFRegMisc::do_nifti_images_match(a_def, NA.get_deformation_field_fwrd(), required_percentage_accuracy))
            throw runtime_error("SIRFRegTransformationAffine::get_as_deformation_field failed.");
        if (!SIRFRegMisc::do_nifti_images_match(b_def, NA.get_deformation_field_fwrd(), required_percentage_accuracy))
            throw runtime_error("SIRFRegTransformationDisplacement::get_as_deformation_field failed.");
        if (!SIRFRegMisc::do_nifti_images_match(c_def, NA.get_deformation_field_fwrd(), required_percentage_accuracy))
            throw runtime_error("SIRFImageDataDeformation::get_as_deformation_field failed.");

        // Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
        mat44 tm_iden = SIRFRegMisc::get_identity_matrix();
        SIRFRegTransformationAffine trans_aff_iden(tm_iden);
        std::vector<std::shared_ptr<SIRFRegTransformation> > vec;
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegTransformationAffine(trans_aff_iden)));
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegTransformationAffine(trans_aff_iden)));
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegTransformationDeformation(c3)));
        SIRFRegTransformationDeformation composed;
        SIRFRegMisc::compose_transformations_into_single_deformation(composed, vec, ref_aladin);
        if (!SIRFRegMisc::do_nifti_images_match(composed.get_as_deformation_field(ref_aladin), NA.get_deformation_field_fwrd(), required_percentage_accuracy))
            throw runtime_error("SIRFRegMisc::compose_transformations_into_single_deformation failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished transformations test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_resample) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting Nifty resample test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        mat44 tm_eye = SIRFRegMisc::get_identity_matrix();
        SIRFRegTransformationAffine tm_iden(tm_eye);
        SIRFRegTransformationAffine tm(NA.get_transformation_matrix_fwrd());
        SIRFRegTransformationDisplacement disp(NA.get_displacement_field_fwrd());
        SIRFRegTransformationDeformation deff(NA.get_deformation_field_fwrd());

        cout << "Testing rigid resample...\n";
        SIRFRegNiftyResample nr1;
        nr1.set_reference_image(ref_aladin);
        nr1.set_floating_image(flo_aladin);
        nr1.set_interpolation_type_to_cubic_spline(); // try different interpolations
        nr1.set_interpolation_type(3); // try different interpolations (cubic)
        nr1.add_transformation_affine(tm_iden);
        nr1.add_transformation_affine(tm);
        nr1.update();
        nr1.get_output().save_to_file(rigid_resample);

        cout << "Testing non-rigid displacement...\n";
        SIRFRegNiftyResample nr2;
        nr2.set_reference_image(ref_aladin);
        nr2.set_floating_image(flo_aladin);
        nr2.set_interpolation_type_to_sinc(); // try different interpolations
        nr2.set_interpolation_type_to_linear(); // try different interpolations
        nr2.add_transformation_disp(disp);
        nr2.update();
        nr2.get_output().save_to_file(nonrigid_resample_disp);

        cout << "Testing non-rigid deformation...\n";
        SIRFRegNiftyResample nr3;
        nr3.set_reference_image(ref_aladin);
        nr3.set_floating_image(flo_aladin);
        nr3.set_interpolation_type_to_nearest_neighbour(); // try different interpolations
        nr3.add_transformation_def(deff);
        nr3.set_interpolation_type_to_linear();
        nr3.update();
        nr3.get_output().save_to_file(nonrigid_resample_def);

        if (!SIRFRegMisc::do_nifti_images_match(NA.get_output(), nr1.get_output(), required_percentage_accuracy))
            throw runtime_error("SIRFRegMisc::compose_transformations_into_single_deformation failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished Nifty resample test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_weighted_mean) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting weighted mean test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        //  Do 3D
        SIRFRegImageWeightedMean3D wm1;
        SIRFImageData im1(stir_nifti);
        SIRFImageData im2(stir_nifti);
        SIRFImageData im3(stir_nifti);
        SIRFImageData im4(stir_nifti);
        im1.fill(1);
        im2.fill(4);
        im3.fill(7);
        im4.fill(6);
        wm1.add_image(im1, 2.F);
        wm1.add_image(im2, 4.F);
        wm1.add_image(im3, 3.F);
        wm1.add_image(im4, 1.F);
        wm1.update();
        wm1.get_output().save_to_file(output_weighted_mean);
        //  Answer should be 4.5, so compare it to that!
        SIRFImageData res(stir_nifti);
        res.fill(4.5F);

        if (!SIRFRegMisc::do_nifti_images_match(wm1.get_output(), res, required_percentage_accuracy))
            throw runtime_error("SIRFRegImageWeightedMean3D failed.");

        //  Do 4D
        SIRFRegImageWeightedMean4D wm2;
        SIRFImageDataDeformation im4D1 = NA.get_deformation_field_fwrd().deep_copy();
        SIRFImageDataDeformation im4D2 = NA.get_deformation_field_fwrd().deep_copy();
        SIRFImageDataDeformation im4D3 = NA.get_deformation_field_fwrd().deep_copy();
        SIRFImageDataDeformation im4D4 = NA.get_deformation_field_fwrd().deep_copy();
        im4D1.fill(1.F);
        im4D2.fill(4.F);
        im4D3.fill(7.F);
        im4D4.fill(6.F);
        wm2.add_image(im4D1, 2.F);
        wm2.add_image(im4D2, 4.F);
        wm2.add_image(im4D3, 3.F);
        wm2.add_image(im4D4, 1.F);
        wm2.update();
        wm2.get_output().save_to_file(output_weighted_mean_def);
        //  Answer should be 4.5, so compare it to that!
        SIRFImageDataDeformation res4D = NA.get_deformation_field_fwrd().deep_copy();
        res4D.fill(4.5);

        if (!SIRFRegMisc::do_nifti_images_match(wm2.get_output(), res4D, required_percentage_accuracy))
            throw runtime_error("SIRFRegImageWeightedMean4D failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished weighted mean test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_stir_to_sirfreg) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting STIR to SIRFReg test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

            // Open stir image
            sirf::PETImageData pet_image_data(stir_nifti);
            SIRFImageData image_data_from_stir(pet_image_data);

            // Compare to nifti IO (if they don't match, you'll see a message but don't throw an error for now)
            SIRFImageData image_data_from_nifti(stir_nifti);
            SIRFRegMisc::do_nifti_images_match(image_data_from_stir, image_data_from_nifti, required_percentage_accuracy);

            // Now fill the stir and sirfreg images with 1 and 100, respectively
            pet_image_data.fill(1.F);
            image_data_from_stir.fill(100.F);

            if (fabs(pet_image_data.data_sptr()->find_max() - image_data_from_stir.get_max()) < 1.e-5F)
                throw runtime_error("STIR & SIRFReg seem to share the same data pointers (their values should be different, but they're the same).");

            // Fill the stir image with the sirfreg
            image_data_from_stir.copy_data_to(pet_image_data);
            if (fabs(pet_image_data.data_sptr()->find_max() - image_data_from_stir.get_max()) > 1.e-5F)
                throw runtime_error("SIRFImageData::copy_data_to failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished STIR to SIRFReg test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    // Error handling
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
