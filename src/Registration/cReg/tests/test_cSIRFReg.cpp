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
#include "NiftiImage3D.h"
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
    string parameter_file_aladin    = examples_path + "/paramFiles/niftyreg_aladin.par";
    string parameter_file_f3d       = examples_path + "/paramFiles/niftyreg_f3d.par";

    // Output filenames
    string save_nifti_image                           = output_path   + "save_NiftiImage";
    string save_nifti_image_3d                        = output_path   + "save_NiftiImage3D";
    string save_nifti_image_3d_tensor_not_split       = output_path   + "save_NiftiImage3DTensor_not_split";
    string save_nifti_image_3d_tensor_split           = output_path   + "save_NiftiImage3DTensor_split_%s";
    string save_nifti_image_3d_deformation_not_split  = output_path   + "save_NiftiImage3DDeformation_not_split";
    string save_nifti_image_3d_deformation_split      = output_path   + "save_NiftiImage3DDeformation_split_%s";
    string save_nifti_image_3d_displacement_not_split = output_path   + "save_NiftiImage3DDisplacement_not_split";
    string save_nifti_image_3d_displacement_split     = output_path   + "save_NiftiImage3DDisplacement_split_%s";
    string aladin_warped            = output_path   + "aladin_warped";
    string f3d_warped               = output_path   + "f3d_warped";
    string TM_fwrd                  = output_path   + "TM_fwrd.txt";
    string TM_back                  = output_path   + "TM_back.txt";
    string aladin_def_fwrd          = output_path   + "aladin_def_fwrd";
    string aladin_def_back          = output_path   + "aladin_def_back_%s";
    string aladin_disp_fwrd         = output_path   + "aladin_disp_fwrd";
    string aladin_disp_back         = output_path   + "aladin_disp_back_%s";
    string f3d_disp_fwrd            = output_path   + "f3d_disp_fwrd";
    string f3d_disp_back            = output_path   + "f3d_disp_back_%s";
    string f3d_def_fwrd             = output_path   + "f3d_def_fwrd";
    string f3d_def_back             = output_path   + "f3d_def_back_%s";
    string rigid_resample           = output_path   + "rigid_resample";
    string nonrigid_resample_disp   = output_path   + "nonrigid_resample_disp";
    string nonrigid_resample_def    = output_path   + "nonrigid_resample_def";
    string output_weighted_mean     = output_path   + "weighted_mean";
    string output_weighted_mean_def = output_path   + "weighted_mean_def";

    string output_stir_nifti        = output_path   + "stir_nifti.nii";

    NiftiImage3D ref_aladin( ref_aladin_filename );
    NiftiImage3D flo_aladin( flo_aladin_filename );
    NiftiImage3D ref_f3d   (   ref_f3d_filename  );
    NiftiImage3D flo_f3d   (   flo_f3d_filename  );

    bool try_misc_functions = true;
    bool try_niftiimage = true;
    bool try_niftiimage3d = true;
    bool try_niftiimage3dtensor = true;
    bool try_niftiimage3ddisplacement = true;
    bool try_niftiimage3ddeformation = true;
    bool try_niftyaladin = true;
    bool try_niftyf3d = true;
    bool try_transformations = true;
    bool try_resample = true;
    bool try_weighted_mean = true;
    bool try_stir_to_sirfreg = true;
    bool try_sirfregmat44 = true;

    if (try_misc_functions) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting misc functions test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // do nifti images match?
        if (ref_aladin != ref_aladin)
            throw runtime_error("Images don't match, but they should.");
        cerr << "\nThe following images intentionally do not match.\n";
        if (ref_aladin == flo_aladin)
            throw runtime_error("Images match, but they shouldn't.");

        // dump from NiftiImage
        ref_aladin.dump_header();
        // dump from multiple images
        NiftiImage::dump_headers({ref_aladin,flo_aladin,ref_f3d,flo_f3d});
        // dump from NiftiImage3DDeformation
        NiftiImage3DDeformation deform;
        deform.create_from_3D_image(ref_aladin);
        deform.dump_header();

        // identity matrix
        SIRFRegMat44 tm_iden = SIRFRegMat44::get_identity();

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished misc functions test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_niftiimage) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting NiftiImage test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // default constructor
        NiftiImage a;

        // Read from file
        NiftiImage b(ref_aladin_filename);

        // Save to file
        b.save_to_file(save_nifti_image);

        // Fill
        b.fill(100);

        // Get max
        if (fabs(b.get_max() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage fill()/get_max() failed.");

        // Get min
        if (fabs(b.get_min() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage fill()/get_min() failed.");

        // Deep copy
        NiftiImage d = b.deep_copy();
        if (&d == &b)
            throw runtime_error("NiftiImage deep_copy failed.");
        if (b != d)
            throw runtime_error("NiftiImage deep_copy failed.");

        // Addition
        NiftiImage3D e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw runtime_error("NiftiImage __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw runtime_error("NiftiImage __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw runtime_error("NiftiImage get_sum() failed.");

        // Add num to image
        NiftiImage3D q = e + 1;
        std::cout << "\nq/e max = " << q.get_max() << "/" << e.get_max() << "\n";
        if (fabs(q.get_max() - (e.get_max() + 1.F)) > 0.0001F)
            throw runtime_error("NiftiImage __add__ val failed.");

        // Subtract num from image
        NiftiImage3D r = e - 1;
        if (fabs(r.get_max() - (e.get_max() - 1.F)) > 0.0001F)
            throw runtime_error("NiftiImage __sub__ val failed.");

        // Multiply image by num
        NiftiImage3D s = e * 10;
        if (fabs(s.get_max() - e.get_max() * 10.F) > 0.0001F)
            throw runtime_error("NiftiImage __mul__ val failed.");

        // Dimensions
        int f[8];
        int g[8] = {3, 64, 64, 64, 1, 1, 1, 1};
        e.get_dimensions(f);
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw runtime_error("NiftiImage get_dimensions() failed.");

        // Test get_element
        int idx[3] = { 1, 2, 3 };
        ref_aladin.get_element(idx);

        // Test get_norm
        if (ref_aladin.get_norm(flo_aladin) < 1.e-7F)
            throw runtime_error("NiftiImage get_norm() failed.");

        // Test get_datatype
        if (strcmp(ref_aladin.get_datatype().c_str(),"NIFTI_TYPE_INT16") != 0)
            throw runtime_error("NiftiImage get_datatype() failed.");

        // Test casting to different datatypes
        vector<string> nifti_types;
        nifti_types.push_back("NIFTI_TYPE_INT16");
        nifti_types.push_back("NIFTI_TYPE_INT32");
        nifti_types.push_back("NIFTI_TYPE_FLOAT32");
        nifti_types.push_back("NIFTI_TYPE_FLOAT64");
        nifti_types.push_back("NIFTI_TYPE_UINT8");
        nifti_types.push_back("NIFTI_TYPE_UINT16");
        nifti_types.push_back("NIFTI_TYPE_UINT32");
        nifti_types.push_back("NIFTI_TYPE_INT64");
        nifti_types.push_back("NIFTI_TYPE_UINT64");
        nifti_types.push_back("NIFTI_TYPE_FLOAT128");
        NiftiImage aa;
        for (int i=0; i<10; ++i) {
            aa = ref_aladin.deep_copy();
            if      (i==0) aa.change_datatype<signed short>();
            else if (i==1) aa.change_datatype<signed int>();
            else if (i==2) aa.change_datatype<float>();
            else if (i==3) aa.change_datatype<double>();
            else if (i==4) aa.change_datatype<unsigned char>();
            else if (i==5) aa.change_datatype<unsigned short>();
            else if (i==6) aa.change_datatype<unsigned int>();
            else if (i==7) aa.change_datatype<signed long long>();
            else if (i==8) aa.change_datatype<unsigned long long>();
            else if (i==9) aa.change_datatype<long double>();

            if (strcmp(aa.get_datatype().c_str(), nifti_types[i].c_str()) != 0)
                throw runtime_error("SIRFRegMisc::change_datatype() failed.");
            if (fabs(aa.get_max() - 255.F) > 1.e-7F)
                throw runtime_error("SIRFRegMisc::change_datatype() failed.");
        }

        // Test dump methods
        q.dump_header();
        NiftiImage::dump_headers({q, s});

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished NiftiImage test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_niftiimage3d) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting NiftiImage3D test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // default constructor
        NiftiImage3D a;

        // Read from file
        NiftiImage3D b(ref_aladin_filename);

        // Save to file
        b.save_to_file(save_nifti_image_3d);

        // Fill
        b.fill(100);

        // Get max
        if (fabs(b.get_max() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage3D fill()/get_max() failed.");

        // Get min
        if (fabs(b.get_min() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage3D fill()/get_min() failed.");

        // Deep copy
        NiftiImage3D d = b.deep_copy();
        if (&d == &b)
            throw runtime_error("NiftiImage3D deep_copy failed.");
        if (d != b)
            throw runtime_error("NiftiImage3D deep_copy failed.");

        // Addition
        NiftiImage3D e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw runtime_error("NiftiImage3D __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw runtime_error("NiftiImage3D __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw runtime_error("NiftiImage3D get_sum() failed.");

        // Dimensions
        int f[8];
        int g[8] = {3, 64, 64, 64, 1, 1, 1, 1};
        e.get_dimensions(f);
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw runtime_error("NiftiImage3D get_dimensions() failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished NiftiImage3D test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_niftiimage3dtensor) {

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting NiftiImage3DTensor test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImage3DTensor from NiftiImage3D
        NiftiImage3DTensor b;
        b.create_from_3D_image(ref_aladin);

        // Save to file
        b.save_to_file(save_nifti_image_3d_tensor_not_split);
        b.save_to_file_split_xyz_components(save_nifti_image_3d_tensor_split);

        // Constructor from file
        NiftiImage3DTensor c(save_nifti_image_3d_tensor_not_split + ".nii");

        // Constructor from single components
        NiftiImage3DTensor h(ref_aladin,ref_aladin,ref_aladin);

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage3DTensor fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage3DTensor fill()/get_min() failed.");

        // Deep copy
        NiftiImage3DTensor d = c.deep_copy();
        if (&d == &c)
            throw runtime_error("NiftiImage3DTensor deep_copy failed.");
        if (d != c)
            throw runtime_error("NiftiImage3DTensor deep_copy failed.");

        // Addition
        NiftiImage3DTensor e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw runtime_error("NiftiImage3DTensor __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw runtime_error("NiftiImage3DTensor __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw runtime_error("NiftiImage3DTensor get_sum() failed.");

        // Dimensions
        int f[8];
        int g[8] = {5, 64, 64, 64, 1, 3, 1, 1};
        e.get_dimensions(f);
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw runtime_error("NiftiImage3DTensor get_dimensions() failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished NiftiImage3DTensor test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_niftiimage3ddisplacement) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting NiftiImage3DDisplacement test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImage3DDisplacement from NiftiImage3D
        NiftiImage3DDisplacement b;
        b.create_from_3D_image(ref_aladin);

        // Save to file
        b.save_to_file(save_nifti_image_3d_displacement_not_split);
        b.save_to_file_split_xyz_components(save_nifti_image_3d_displacement_split);

        // Constructor from file
        NiftiImage3DDisplacement c(save_nifti_image_3d_displacement_not_split + ".nii");

        // Constructor from tensor
        NiftiImage3DTensor x(save_nifti_image_3d_displacement_not_split + ".nii");
        NiftiImage3DDisplacement y(x);

        // Constructor from general
        NiftiImage q(save_nifti_image_3d_displacement_not_split + ".nii");
        NiftiImage3DDisplacement r(q);

        // Constructor from single components
        NiftiImage3DDisplacement h(ref_aladin,ref_aladin,ref_aladin);

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage3DDisplacement fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage3DDisplacement fill()/get_min() failed.");

        // Deep copy
        NiftiImage3DDisplacement d = c.deep_copy();
        if (&d == &c)
            throw runtime_error("NiftiImage3DDisplacement deep_copy failed.");
        if (d != c)
            throw runtime_error("NiftiImage3DDisplacement deep_copy failed.");

        // Addition
        NiftiImage3DDisplacement e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw runtime_error("NiftiImage3DDisplacement __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw runtime_error("NiftiImage3DDisplacement __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw runtime_error("NiftiImage3DDisplacement get_sum() failed.");

        // Dimensions
        int f[8];
        int g[8] = {5, 64, 64, 64, 1, 3, 1, 1};
        e.get_dimensions(f);
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw runtime_error("NiftiImage3DDisplacement get_dimensions() failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished NiftiImage3DDisplacement test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_niftiimage3ddeformation) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting NiftiImage3DDeformation test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImage3DDeformation from NiftiImage3D
        NiftiImage3DDeformation b;
        b.create_from_3D_image(ref_aladin);

        // Save to file
        b.save_to_file(save_nifti_image_3d_deformation_not_split);
        b.save_to_file_split_xyz_components(save_nifti_image_3d_deformation_split);

        // Constructor from file
        NiftiImage3DDeformation c(save_nifti_image_3d_deformation_not_split + ".nii");

        // Constructor from tensor
        NiftiImage3DTensor x(save_nifti_image_3d_deformation_not_split + ".nii");
        NiftiImage3DDeformation y(x);

        // Constructor from general
        NiftiImage q(save_nifti_image_3d_deformation_not_split + ".nii");
        NiftiImage3DDeformation r(q);

        // Constructor from single components
        NiftiImage3DDeformation h(ref_aladin,ref_aladin,ref_aladin);

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage3DDeformation fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw runtime_error("NiftiImage3DDeformation fill()/get_min() failed.");

        // Deep copy
        NiftiImage3DDeformation d = c.deep_copy();
        if (&d == &c)
            throw runtime_error("NiftiImage3DDeformation deep_copy failed.");
        if (d != c)
            throw runtime_error("NiftiImage3DDeformation deep_copy failed.");

        // Addition
        NiftiImage3DDeformation e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw runtime_error("NiftiImage3DDeformation __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw runtime_error("NiftiImage3DDeformation __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw runtime_error("NiftiImage3DDeformation get_sum() failed.");

        // Dimensions
        int f[8];
        int g[8] = {5, 64, 64, 64, 1, 3, 1, 1};
        e.get_dimensions(f);
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw runtime_error("NiftiImage3DDeformation get_dimensions() failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished NiftiImage3DDeformation test.\n";
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
        NA.get_transformation_matrix_fwrd().save_to_file(       TM_fwrd      );
        NA.get_transformation_matrix_back().save_to_file(       TM_back      );
        NA.get_displacement_field_fwrd().save_to_file(aladin_disp_fwrd);
        NA.get_displacement_field_back().save_to_file_split_xyz_components(aladin_disp_back);
        NA.get_deformation_field_fwrd().save_to_file(aladin_def_fwrd);
        NA.get_deformation_field_back().save_to_file_split_xyz_components(aladin_def_back);

        // Get outputs
        NiftiImage3D warped = NA.get_output();
        NiftiImage3DTensor def_fwrd  = NA.get_deformation_field_fwrd();
        NiftiImage3DTensor def_back  = NA.get_deformation_field_back();
        NiftiImage3DTensor disp_fwrd = NA.get_displacement_field_fwrd();
        NiftiImage3DTensor disp_back = NA.get_displacement_field_back();

        // Fwrd TM
        SIRFRegMat44 fwrd_tm = NA.get_transformation_matrix_fwrd();
        fwrd_tm.print();

        // Back TM
        SIRFRegMat44 back_tm = NA.get_transformation_matrix_back();
        back_tm.print();

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
        NF.get_deformation_field_fwrd().save_to_file (f3d_def_fwrd);
        NF.get_deformation_field_back().save_to_file_split_xyz_components(f3d_def_back);
        NF.get_displacement_field_fwrd().save_to_file(f3d_disp_fwrd);
        NF.get_displacement_field_back().save_to_file_split_xyz_components(f3d_disp_back);

        // Get outputs
        NiftiImage3D warped = NF.get_output();
        NiftiImage3DTensor def_fwrd  = NF.get_deformation_field_fwrd();
        NiftiImage3DTensor def_back  = NF.get_deformation_field_back();
        NiftiImage3DTensor disp_fwrd = NF.get_displacement_field_fwrd();
        NiftiImage3DTensor disp_back = NF.get_displacement_field_back();

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished Nifty f3d test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_transformations) {
        if (!try_niftyaladin)
            throw runtime_error("This test requires you to have run aladin.");
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting transformations test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        // Affine
        cout << "\nTesting affine...\n";
        SIRFRegMat44 a1;
        SIRFRegMat44 a2(TM_fwrd);
        SIRFRegMat44 a3(NA.get_transformation_matrix_fwrd());

        // Displacement
        cout << "\nTesting displacement...\n";
        NiftiImage3DDisplacement b3(NA.get_displacement_field_fwrd());

        // Deformation
        cout << "\nTesting deformation...\n";
        NiftiImage3DDeformation c3(NA.get_deformation_field_fwrd());

        // Get as deformations
        NiftiImage3DDeformation a_def = a3.get_as_deformation_field(ref_aladin);
        NiftiImage3DDeformation b_def = b3.get_as_deformation_field(ref_aladin);
        NiftiImage3DDeformation c_def = c3.get_as_deformation_field(ref_aladin);
        if (a_def != NA.get_deformation_field_fwrd())
            throw runtime_error("SIRFRegMat44::get_as_deformation_field failed.");
        if (b_def != NA.get_deformation_field_fwrd())
            throw runtime_error("NiftiImage3DDisplacement::get_as_deformation_field failed.");
        if (c_def != NA.get_deformation_field_fwrd())
            throw runtime_error("NiftiImage3DDeformation::get_as_deformation_field failed.");

        // Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
        SIRFRegMat44 tm_iden = SIRFRegMat44::get_identity();
        SIRFRegMat44 trans_aff_iden(tm_iden);
        std::vector<std::shared_ptr<SIRFRegTransformation> > vec;
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegMat44(trans_aff_iden)));
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegMat44(trans_aff_iden)));
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new NiftiImage3DDeformation(c3)));
        NiftiImage3DDeformation composed =
                NiftiImage3DDeformation::compose_single_deformation(vec, ref_aladin);
        if (composed.get_as_deformation_field(ref_aladin) != NA.get_deformation_field_fwrd())
            throw runtime_error("SIRFRegMisc::compose_transformations_into_single_deformation failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished transformations test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_resample) {
        if (!try_niftyaladin)
            throw runtime_error("This test requires you to have run aladin.");
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting Nifty resample test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        SIRFRegMat44 tm_eye = SIRFRegMat44::get_identity();
        SIRFRegMat44 tm_iden(tm_eye);
        SIRFRegMat44 tm(NA.get_transformation_matrix_fwrd());
        NiftiImage3DDisplacement disp(NA.get_displacement_field_fwrd());
        NiftiImage3DDeformation deff(NA.get_deformation_field_fwrd());

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

        if (NA.get_output() != nr1.get_output())
            throw runtime_error("SIRFRegMisc::compose_transformations_into_single_deformation failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished Nifty resample test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_weighted_mean) {
        if (!try_niftyaladin)
            throw runtime_error("This test requires you to have run aladin.");
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting weighted mean test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

        //  Do 3D
        SIRFRegImageWeightedMean wm1;
        // Change to float to avoid rounding errors
        NiftiImage3D ref_aladin_float(ref_aladin);
        ref_aladin_float.change_datatype<float>();
        NiftiImage3D im1 = ref_aladin_float.deep_copy();
        NiftiImage3D im2 = ref_aladin_float.deep_copy();
        NiftiImage3D im3 = ref_aladin_float.deep_copy();
        NiftiImage3D im4 = ref_aladin_float.deep_copy();
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
        NiftiImage3D res = ref_aladin_float;
        res.fill(4.5F);

        if (wm1.get_output() != res)
            throw runtime_error("SIRFRegImageWeightedMean3D failed.");

        //  Do 4D
        SIRFRegImageWeightedMean wm2;
        NiftiImage3DTensor im4D1 = NA.get_deformation_field_fwrd().deep_copy();
        NiftiImage3DTensor im4D2 = NA.get_deformation_field_fwrd().deep_copy();
        NiftiImage3DTensor im4D3 = NA.get_deformation_field_fwrd().deep_copy();
        NiftiImage3DTensor im4D4 = NA.get_deformation_field_fwrd().deep_copy();
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
        NiftiImage3DTensor res4D = NA.get_deformation_field_fwrd().deep_copy();
        res4D.fill(4.5);

        if (wm2.get_output() != res4D)
            throw runtime_error("SIRFRegImageWeightedMean3DTensor failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished weighted mean test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_stir_to_sirfreg) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting STIR to SIRFReg test...\n";
        cout << "//------------------------------------------------------------------------ //\n";

            // Open stir image
            sirf::PETImageData pet_image_data(ref_aladin_filename);
            NiftiImage3D image_data_from_stir(pet_image_data);

            // Now fill the stir and sirfreg images with 1 and 100, respectively
            pet_image_data.fill(1.F);
            image_data_from_stir.fill(100.F);

            if (fabs(pet_image_data.data_sptr()->find_max() - image_data_from_stir.get_max()) < 1.e-5F)
                throw runtime_error("STIR & SIRFReg seem to share the same data pointers (their values should be different, but they're the same).");

            // Fill the stir image with the sirfreg
            image_data_from_stir.copy_data_to(pet_image_data);
            if (fabs(pet_image_data.data_sptr()->find_max() - image_data_from_stir.get_max()) > 1.e-5F)
                throw runtime_error("NiftiImage3D::copy_data_to failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished STIR to SIRFReg test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    if (try_sirfregmat44) {
        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Starting SIRFRegMat44 test...\n";
        cout << "//------------------------------------------------------------------------ //\n";
        if (!try_niftyaladin)
            throw runtime_error("This test requires you to have run aladin.");

        // Construct from file
        SIRFRegMat44 a(TM_fwrd);

        // Multiply fwrd and inverse, should equal identity
        SIRFRegMat44 b = NA.get_transformation_matrix_fwrd();
        SIRFRegMat44 c = NA.get_transformation_matrix_back();
        SIRFRegMat44 d = b * c;
        SIRFRegMat44 e = SIRFRegMat44::get_identity();
        if (d != e)
            throw runtime_error("SIRFRegMat44::mult/comparison failed.");

        d.fill(3.F);
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                if ((d[i][j] - 3.F) > 1.e-7F)
                    throw runtime_error("SIRFRegMat44::fill/operator[] failed.");

        if (d.get_determinant() > 1.e-7F)
            throw runtime_error("SIRFRegMat44::get_determinant failed.");
        if (e.get_determinant() - 1.F > 1.e-7F)
            throw runtime_error("SIRFRegMat44::get_determinant failed.");

        cout << "// ----------------------------------------------------------------------- //\n";
        cout << "//                  Finished SIRFRegMat44 test.\n";
        cout << "//------------------------------------------------------------------------ //\n";
    }

    // Error handling
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
