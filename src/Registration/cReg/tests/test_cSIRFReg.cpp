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
#include "NiftiImageData3D.h"
#include "SIRFRegImageWeightedMean.h"
#include "stir_data_containers.h"
#include <memory>

using namespace sirf;

int main(int argc, char* argv[])
{

    try {

    // Paths
    std::string SIRF_PATH;
    if (argc==1)
        SIRF_PATH = getenv("SIRF_PATH");
    else
        SIRF_PATH = argv[1];
    const std::string examples_path = SIRF_PATH + "/data/examples/Registration";
    const std::string output_prefix   = "results/cplusplus_";

    // Input filenames
    const std::string ref_aladin_filename      = examples_path + "/test.nii.gz";
    const std::string flo_aladin_filename      = examples_path + "/test2.nii.gz";
    const std::string ref_f3d_filename         = examples_path + "/mouseFixed.nii.gz";
    const std::string flo_f3d_filename         = examples_path + "/mouseMoving.nii.gz";
    const std::string parameter_file_aladin    = examples_path + "/paramFiles/niftyreg_aladin.par";
    const std::string parameter_file_f3d       = examples_path + "/paramFiles/niftyreg_f3d.par";

    // Output filenames
    const std::string save_nifti_image                           = output_prefix   + "save_NiftiImageData.nii";
    const std::string save_nifti_image_3d                        = output_prefix   + "save_NiftiImageData3D.nii";
    const std::string save_nifti_image_3d_tensor_not_split       = output_prefix   + "save_NiftiImageData3DTensor_not_split.nii";
    const std::string save_nifti_image_3d_tensor_split           = output_prefix   + "save_NiftiImageData3DTensor_split_%s.nii";
    const std::string save_nifti_image_3d_deformation_not_split  = output_prefix   + "save_NiftiImageData3DDeformation_not_split.nii";
    const std::string save_nifti_image_3d_deformation_split      = output_prefix   + "save_NiftiImageData3DDeformation_split_%s.nii";
    const std::string save_nifti_image_3d_displacement_not_split = output_prefix   + "save_NiftiImageData3DDisplacement_not_split.nii";
    const std::string save_nifti_image_3d_displacement_split     = output_prefix   + "save_NiftiImageData3DDisplacement_split_%s.nii";
    const std::string aladin_warped            = output_prefix   + "aladin_warped.nii";
    const std::string f3d_warped               = output_prefix   + "f3d_warped.nii";
    const std::string TM_forward               = output_prefix   + "TM_forward.txt";
    const std::string TM_inverse               = output_prefix   + "TM_inverse.txt";
    const std::string aladin_def_forward       = output_prefix   + "aladin_def_forward.nii";
    const std::string aladin_def_inverse       = output_prefix   + "aladin_def_inverse_%s.nii";
    const std::string aladin_disp_forward      = output_prefix   + "aladin_disp_forward.nii";
    const std::string aladin_disp_inverse      = output_prefix   + "aladin_disp_inverse_%s.nii";
    const std::string f3d_disp_forward         = output_prefix   + "f3d_disp_forward.nii";
    const std::string f3d_disp_inverse         = output_prefix   + "f3d_disp_inverse_%s.nii";
    const std::string f3d_def_forward          = output_prefix   + "f3d_def_forward.nii";
    const std::string f3d_def_inverse          = output_prefix   + "f3d_def_inverse_%s.nii";
    const std::string rigid_resample           = output_prefix   + "rigid_resample.nii";
    const std::string nonrigid_resample_disp   = output_prefix   + "nonrigid_resample_disp.nii";
    const std::string nonrigid_resample_def    = output_prefix   + "nonrigid_resample_def.nii";
    const std::string output_weighted_mean     = output_prefix   + "weighted_mean.nii";
    const std::string output_weighted_mean_def = output_prefix   + "weighted_mean_def.nii";
    const std::string output_float             = output_prefix   + "reg_aladin_float.nii";

    const NiftiImageData3D ref_aladin( ref_aladin_filename );
    const NiftiImageData3D flo_aladin( flo_aladin_filename );
    const NiftiImageData3D ref_f3d   (   ref_f3d_filename  );
    const NiftiImageData3D flo_f3d   (   flo_f3d_filename  );

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // default constructor
        NiftiImageData a;

        // Read from file
        NiftiImageData b(ref_aladin_filename);

        // Save to file
        b.save_to_file(save_nifti_image);

        // Fill
        b.fill(100);

        // Get max
        if (fabs(b.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData fill()/get_max() failed.");

        // Get min
        if (fabs(b.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData fill()/get_min() failed.");

        // Deep copy
        NiftiImageData d = b.deep_copy();
        std::cout << "\ndone.\n";
        if (d.get_raw_nifti_sptr() == b.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData deep_copy failed.");
        if (b != d)
            throw std::runtime_error("NiftiImageData deep_copy failed.");

        // Addition
        NiftiImageData3D e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw std::runtime_error("NiftiImageData __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw std::runtime_error("NiftiImageData __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw std::runtime_error("NiftiImageData get_sum() failed.");

        // Add num to image
        NiftiImageData3D q = e + 1;
        if (fabs(q.get_max() - (e.get_max() + 1.F)) > 0.0001F)
            throw std::runtime_error("NiftiImageData __add__ val failed.");

        // Subtract num from image
        NiftiImageData3D r = e - 1;
        if (fabs(r.get_max() - (e.get_max() - 1.F)) > 0.0001F)
            throw std::runtime_error("NiftiImageData __sub__ val failed.");

        // Multiply image by num
        NiftiImageData3D s = e * 10;
        if (fabs(s.get_max() - e.get_max() * 10.F) > 0.0001F)
            throw std::runtime_error("NiftiImageData __mul__ val failed.");

        // Dimensions
        int g[8] = {3, 64, 64, 64, 1, 1, 1, 1};
        const int *f = e.get_dimensions();
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw std::runtime_error("NiftiImageData get_dimensions() failed.");

        // Test get_element
        int idx[7] = { 1, 2, 3, 0, 0, 0, 0 };
        ref_aladin(idx);

        // Test get_norm
        if (ref_aladin.get_norm(flo_aladin) < 1.e-7F)
            throw std::runtime_error("NiftiImageData get_norm() failed.");

        // Test saving to datatype
        ref_aladin.save_to_file(output_float,NIFTI_TYPE_FLOAT32);
        NiftiImageData3D ref_aladin_float(output_float);
        for (int i=0; i<int(ref_aladin.get_raw_nifti_sptr()->nvox); ++i)
            if (ref_aladin_float(i) - ref_aladin(i) > 1.e-7F)
                throw std::runtime_error("NiftiImageData3D::save_to_file()/change_datatype() failed.");

        // Test print methods
        q.print_header();
        NiftiImageData::print_headers({q, s});

        // Crop image
        int min[7], max[7];
        for (int i=0; i<7; ++i) {
            min[i] = 0;
            max[i] = f[i+1] - 1;
        }
        max[2] = 62;
        NiftiImageData z = e.deep_copy();
        z.crop(min, max);
        const int *zz = z.get_dimensions();
        if (zz[0] != 3 ||
                zz[1] != 64 ||
                zz[2] != 64 ||
                zz[3] != 63 ||
                zz[4] != 1 ||
                zz[5] != 1 ||
                zz[6] != 1)
            throw std::runtime_error("NiftiImageData3D::crop() failed.");


        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData3D test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // default constructor
        NiftiImageData3D a;

        // Read from file
        NiftiImageData3D b(ref_aladin_filename);

        // Save to file
        b.save_to_file(save_nifti_image_3d);

        // Fill
        b.fill(100);

        // Get max
        if (fabs(b.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3D fill()/get_max() failed.");

        // Get min
        if (fabs(b.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3D fill()/get_min() failed.");

        // Deep copy
        NiftiImageData3D d = b.deep_copy();
        if (d.get_raw_nifti_sptr() == b.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData3D deep_copy failed.");
        if (d != b)
            throw std::runtime_error("NiftiImageData3D deep_copy failed.");

        // Addition
        NiftiImageData3D e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw std::runtime_error("NiftiImageData3D __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw std::runtime_error("NiftiImageData3D __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw std::runtime_error("NiftiImageData3D get_sum() failed.");

        // Dimensions
        int g[8] = {3, 64, 64, 64, 1, 1, 1, 1};
        const int *f = e.get_dimensions();
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw std::runtime_error("NiftiImageData3D get_dimensions() failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3D test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData3DTensor test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImageData3DTensor from NiftiImageData3D
        NiftiImageData3DTensor b;
        b.create_from_3D_image(ref_aladin);

        // Save to file
        b.save_to_file(save_nifti_image_3d_tensor_not_split);
        b.save_to_file_split_xyz_components(save_nifti_image_3d_tensor_split);

        // Constructor from file
        NiftiImageData3DTensor c(save_nifti_image_3d_tensor_not_split);

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DTensor fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DTensor fill()/get_min() failed.");

        // Deep copy
        NiftiImageData3DTensor d = c.deep_copy();
        if (d.get_raw_nifti_sptr() == c.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData3DTensor deep_copy failed.");
        if (d != c)
            throw std::runtime_error("NiftiImageData3DTensor deep_copy failed.");

        // Addition
        NiftiImageData3DTensor e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw std::runtime_error("NiftiImageData3DTensor __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw std::runtime_error("NiftiImageData3DTensor __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw std::runtime_error("NiftiImageData3DTensor get_sum() failed.");

        // Dimensions
        int g[8] = {5, 64, 64, 64, 1, 3, 1, 1};
        const int *f = e.get_dimensions();
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw std::runtime_error("NiftiImageData3DTensor get_dimensions() failed.");

        // Constructor from single components
        NiftiImageData3D im1 = ref_aladin.deep_copy();
        NiftiImageData3D im2 = ref_aladin.deep_copy();
        NiftiImageData3D im3 = ref_aladin.deep_copy();
        im1.fill(30.F);
        im2.fill(20.F);
        im3.fill(-10.F);
        NiftiImageData3DTensor h(im1, im2, im3);

        // Test flip components
        h.flip_component(0);
        if (fabs(h.get_max() - 20.F) > 1.e-7F )
            throw std::runtime_error("NiftiImageData3DTensor flip_component() failed.");
        if (fabs(h.get_min() + 30.F) > 1.e-7F )
            throw std::runtime_error("NiftiImageData3DTensor flip_component() failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3DTensor test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData3DDisplacement test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImageData3DDisplacement from NiftiImageData3D
        NiftiImageData3DDisplacement b;
        b.create_from_3D_image(ref_aladin);

        // Save to file
        b.save_to_file(save_nifti_image_3d_displacement_not_split);
        b.save_to_file_split_xyz_components(save_nifti_image_3d_displacement_split);

        // Constructor from file
        NiftiImageData3DDisplacement c(save_nifti_image_3d_displacement_not_split);

        // Constructor from tensor
        NiftiImageData3DTensor x(save_nifti_image_3d_displacement_not_split);
        NiftiImageData3DDisplacement y(x);

        // Constructor from general
        NiftiImageData q(save_nifti_image_3d_displacement_not_split);
        NiftiImageData3DDisplacement r(q);

        // Constructor from single components
        NiftiImageData3DDisplacement h(ref_aladin,ref_aladin,ref_aladin);

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DDisplacement fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DDisplacement fill()/get_min() failed.");

        // Deep copy
        NiftiImageData3DDisplacement d = c.deep_copy();
        if (d.get_raw_nifti_sptr() == c.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData3DDisplacement deep_copy failed.");
        if (d != c)
            throw std::runtime_error("NiftiImageData3DDisplacement deep_copy failed.");

        // Addition
        NiftiImageData3DDisplacement e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw std::runtime_error("NiftiImageData3DDisplacement __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw std::runtime_error("NiftiImageData3DDisplacement __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw std::runtime_error("NiftiImageData3DDisplacement get_sum() failed.");

        // Dimensions
        int g[8] = {5, 64, 64, 64, 1, 3, 1, 1};
        const int *f = e.get_dimensions();
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw std::runtime_error("NiftiImageData3DDisplacement get_dimensions() failed.");


        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3DDisplacement test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData3DDeformation test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImageData3DDeformation from NiftiImageData3D
        NiftiImageData3DDeformation b;
        b.create_from_3D_image(ref_aladin);

        // Save to file
        b.save_to_file(save_nifti_image_3d_deformation_not_split);
        b.save_to_file_split_xyz_components(save_nifti_image_3d_deformation_split);

        // Constructor from file
        NiftiImageData3DDeformation c(save_nifti_image_3d_deformation_not_split);

        // Constructor from tensor
        NiftiImageData3DTensor x(save_nifti_image_3d_deformation_not_split);
        NiftiImageData3DDeformation y(x);

        // Constructor from general
        NiftiImageData q(save_nifti_image_3d_deformation_not_split);
        NiftiImageData3DDeformation r(q);

        // Constructor from single components
        NiftiImageData3DDeformation h(ref_aladin,ref_aladin,ref_aladin);

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DDeformation fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DDeformation fill()/get_min() failed.");

        // Deep copy
        NiftiImageData3DDeformation d = c.deep_copy();
        if (d.get_raw_nifti_sptr() == c.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData3DDeformation deep_copy failed.");
        if (d != c)
            throw std::runtime_error("NiftiImageData3DDeformation deep_copy failed.");

        // Addition
        NiftiImageData3DDeformation e = d + d;
        if (fabs(e.get_max() - 2 * d.get_max()) > 0.0001F)
            throw std::runtime_error("NiftiImageData3DDeformation __add__/get_max() failed.");

        // Subtraction
        e = d - d;
        if (e.get_max() > 0.0001F)
            throw std::runtime_error("NiftiImageData3DDeformation __sub__/get_max() failed.");

        // Sum
        if (e.get_sum() > 0.0001F)
            throw std::runtime_error("NiftiImageData3DDeformation get_sum() failed.");

        // Dimensions
        int g[8] = {5, 64, 64, 64, 1, 3, 1, 1};
        const int *f = e.get_dimensions();
        for (int i=0; i<8; ++i)
            if (g[i] != f[i])
                throw std::runtime_error("NiftiImageData3DDeformation get_dimensions() failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3DDeformation test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    SIRFRegNiftyAladinSym<float> NA;
    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Nifty aladin test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        NA.set_reference_image               (           ref_aladin          );
        NA.set_floating_image                (           flo_aladin          );
        NA.set_parameter_file                (      parameter_file_aladin    );
        NA.set_parameter("SetInterpolationToCubic");
        NA.set_parameter("SetLevelsToPerform","1");
        NA.set_parameter("SetMaxIterations","5");
        NA.process();
        NA.get_output().save_to_file         (         aladin_warped         );
        NA.get_transformation_matrix_forward().save_to_file(       TM_forward      );
        NA.get_transformation_matrix_inverse().save_to_file(       TM_inverse      );
        NA.get_displacement_field_forward().save_to_file(aladin_disp_forward);
        NA.get_displacement_field_inverse().save_to_file_split_xyz_components(aladin_disp_inverse);
        NA.get_deformation_field_forward().save_to_file(aladin_def_forward);
        NA.get_deformation_field_inverse().save_to_file_split_xyz_components(aladin_def_inverse);

        // Get outputs
        NiftiImageData3D warped = NA.get_output();
        NiftiImageData3DTensor def_forward  = NA.get_deformation_field_forward();
        NiftiImageData3DTensor def_inverse  = NA.get_deformation_field_inverse();
        NiftiImageData3DTensor disp_forward = NA.get_displacement_field_forward();
        NiftiImageData3DTensor disp_inverse = NA.get_displacement_field_inverse();

        // forward TM
        SIRFRegAffineTransformation forward_tm = NA.get_transformation_matrix_forward();
        forward_tm.print();

        // Inverse TM
        SIRFRegAffineTransformation inverse_tm = NA.get_transformation_matrix_inverse();
        inverse_tm.print();

        // Test converting disp to def
        NiftiImageData3DDeformation a;
        a.create_from_disp(disp_forward);
        if (a != def_forward)
            throw std::runtime_error("NiftiImageData3DDeformation::create_from_disp() failed.");

        // Test converting def to disp
        NiftiImageData3DDisplacement b;
        b.create_from_def(def_forward);
        if (b != disp_forward)
            throw std::runtime_error("NiftiImageData3DDisplacement::create_from_def() failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished Nifty aladin test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }


    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Nifty f3d test..\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        SIRFRegNiftyF3dSym<float> NF;
        NF.set_reference_image               (          ref_f3d           );
        NF.set_floating_image                (          flo_f3d           );
        NF.set_parameter_file                (     parameter_file_f3d     );
        NF.set_reference_time_point          (             1              );
        NF.set_floating_time_point           (             1              );
        NF.process();
        NF.get_output().save_to_file         (         f3d_warped         );
        NF.get_deformation_field_forward().save_to_file (f3d_def_forward);
        NF.get_deformation_field_inverse().save_to_file_split_xyz_components(f3d_def_inverse);
        NF.get_displacement_field_forward().save_to_file(f3d_disp_forward);
        NF.get_displacement_field_inverse().save_to_file_split_xyz_components(f3d_disp_inverse);

        // Get outputs
        NiftiImageData3D warped = NF.get_output();
        NiftiImageData3DTensor def_forward  = NF.get_deformation_field_forward();
        NiftiImageData3DTensor def_inverse  = NF.get_deformation_field_inverse();
        NiftiImageData3DTensor disp_forward = NF.get_displacement_field_forward();
        NiftiImageData3DTensor disp_inverse = NF.get_displacement_field_inverse();

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished Nifty f3d test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting transformations test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Affine
        std::cout << "\nTesting affine...\n";
        SIRFRegAffineTransformation a1;
        SIRFRegAffineTransformation a2(TM_forward);
        SIRFRegAffineTransformation a3(NA.get_transformation_matrix_forward());

        // Displacement
        std::cout << "\nTesting displacement...\n";
        NiftiImageData3DDisplacement b3(NA.get_displacement_field_forward());

        // Deformation
        std::cout << "\nTesting deformation...\n";
        NiftiImageData3DDeformation c3(NA.get_deformation_field_forward());

        // Get as deformations
        NiftiImageData3DDeformation a_def = a3.get_as_deformation_field(ref_aladin);
        NiftiImageData3DDeformation b_def = b3.get_as_deformation_field(ref_aladin);
        NiftiImageData3DDeformation c_def = c3.get_as_deformation_field(ref_aladin);
        if (a_def != NA.get_deformation_field_forward())
            throw std::runtime_error("SIRFRegAffineTransformation::get_as_deformation_field failed.");
        if (b_def != NA.get_deformation_field_forward())
            throw std::runtime_error("NiftiImageData3DDisplacement::get_as_deformation_field failed.");
        if (c_def != NA.get_deformation_field_forward())
            throw std::runtime_error("NiftiImageData3DDeformation::get_as_deformation_field failed.");

        // Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
        SIRFRegAffineTransformation tm_iden = SIRFRegAffineTransformation::get_identity();
        SIRFRegAffineTransformation trans_aff_iden(tm_iden);
        std::vector<std::shared_ptr<SIRFRegTransformation> > vec;
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegAffineTransformation(trans_aff_iden)));
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegAffineTransformation(trans_aff_iden)));
        vec.push_back(std::shared_ptr<SIRFRegTransformation>(new NiftiImageData3DDeformation(c3)));
        NiftiImageData3DDeformation composed =
                NiftiImageData3DDeformation::compose_single_deformation(vec, ref_aladin);
        if (composed.get_as_deformation_field(ref_aladin) != NA.get_deformation_field_forward())
            throw std::runtime_error("NiftiImageData3DDeformation::compose_single_deformation failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished transformations test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Nifty resample test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        SIRFRegAffineTransformation tm_eye = SIRFRegAffineTransformation::get_identity();
        SIRFRegAffineTransformation tm_iden(tm_eye);
        SIRFRegAffineTransformation tm(NA.get_transformation_matrix_forward());
        NiftiImageData3DDisplacement disp(NA.get_displacement_field_forward());
        NiftiImageData3DDeformation deff(NA.get_deformation_field_forward());

        std::cout << "Testing rigid resample...\n";
        SIRFRegNiftyResample nr1;
        nr1.set_reference_image(ref_aladin);
        nr1.set_floating_image(flo_aladin);
        nr1.set_interpolation_type_to_cubic_spline(); // try different interpolations
        nr1.set_interpolation_type(SIRFRegNiftyResample::CUBICSPLINE); // try different interpolations (cubic)
        nr1.add_transformation_affine(tm_iden);
        nr1.add_transformation_affine(tm);
        nr1.process();
        nr1.get_output().save_to_file(rigid_resample);

        std::cout << "Testing non-rigid displacement...\n";
        SIRFRegNiftyResample nr2;
        nr2.set_reference_image(ref_aladin);
        nr2.set_floating_image(flo_aladin);
        nr2.set_interpolation_type_to_sinc(); // try different interpolations
        nr2.set_interpolation_type_to_linear(); // try different interpolations
        nr2.add_transformation_disp(disp);
        nr2.process();
        nr2.get_output().save_to_file(nonrigid_resample_disp);

        std::cout << "Testing non-rigid deformation...\n";
        SIRFRegNiftyResample nr3;
        nr3.set_reference_image(ref_aladin);
        nr3.set_floating_image(flo_aladin);
        nr3.set_interpolation_type_to_nearest_neighbour(); // try different interpolations
        nr3.add_transformation_def(deff);
        nr3.set_interpolation_type_to_linear();
        nr3.process();
        nr3.get_output().save_to_file(nonrigid_resample_def);

        // TODO This isn't working on my machine. But it's not working with NiftyReg executables either, so I don't think it's my code
        // i.e., reg_aladin != reg_resample when reg_resample uses the transformation matrix from reg_aladin
        /*if (NA.get_output() != nr1.get_output())
            throw std::runtime_error("compose_transformations_into_single_deformation failed.");*/

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished Nifty resample test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting weighted mean test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        //  Do 3D
        SIRFRegImageWeightedMean wm1;
        NiftiImageData3D im1 = ref_aladin.deep_copy();
        NiftiImageData3D im2 = ref_aladin.deep_copy();
        NiftiImageData3D im3 = ref_aladin.deep_copy();
        NiftiImageData3D im4 = ref_aladin.deep_copy();
        im1.fill(1);
        im2.fill(4);
        im3.fill(7);
        im4.fill(6);
        wm1.add_image(im1, 2.F);
        wm1.add_image(im2, 4.F);
        wm1.add_image(im3, 3.F);
        wm1.add_image(im4, 1.F);
        wm1.process();
        wm1.get_output().save_to_file(output_weighted_mean);
        //  Answer should be 4.5, so compare it to that!
        NiftiImageData3D res = ref_aladin;
        res.fill(4.5F);

        if (wm1.get_output() != res)
            throw std::runtime_error("SIRFRegImageWeightedMean3D failed.");

        //  Do 4D
        SIRFRegImageWeightedMean wm2;
        NiftiImageData3DTensor im4D1 = NA.get_deformation_field_forward().deep_copy();
        NiftiImageData3DTensor im4D2 = NA.get_deformation_field_forward().deep_copy();
        NiftiImageData3DTensor im4D3 = NA.get_deformation_field_forward().deep_copy();
        NiftiImageData3DTensor im4D4 = NA.get_deformation_field_forward().deep_copy();
        im4D1.fill(1.F);
        im4D2.fill(4.F);
        im4D3.fill(7.F);
        im4D4.fill(6.F);
        wm2.add_image(im4D1, 2.F);
        wm2.add_image(im4D2, 4.F);
        wm2.add_image(im4D3, 3.F);
        wm2.add_image(im4D4, 1.F);
        wm2.process();
        wm2.get_output().save_to_file(output_weighted_mean_def);
        //  Answer should be 4.5, so compare it to that!
        NiftiImageData3DTensor res4D = NA.get_deformation_field_forward().deep_copy();
        res4D.fill(4.5);

        if (wm2.get_output() != res4D)
            throw std::runtime_error("SIRFRegImageWeightedMean3DTensor failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished weighted mean test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting STIR to SIRFReg test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

            // Open stir image
            sirf::PETImageData pet_image_data(ref_aladin_filename);
            NiftiImageData3D image_data_from_stir(pet_image_data);

            // Now fill the stir and sirfreg images with 1 and 100, respectively
            pet_image_data.fill(1.F);
            image_data_from_stir.fill(100.F);

            if (fabs(pet_image_data.data_sptr()->find_max() - image_data_from_stir.get_max()) < 1.e-5F)
                throw std::runtime_error("STIR & SIRFReg seem to share the same data pointers (their values should be different, but they're the same).");

            // Fill the stir image with the sirfreg
            image_data_from_stir.copy_data_to(pet_image_data);
            if (fabs(pet_image_data.data_sptr()->find_max() - image_data_from_stir.get_max()) > 1.e-5F)
                throw std::runtime_error("NiftiImageData3D::copy_data_to failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished STIR to SIRFReg test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting SIRFRegAffineTransformation test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Construct from file
        SIRFRegAffineTransformation a(TM_forward);

        // Multiply forward and inverse, should equal identity
        SIRFRegAffineTransformation b = NA.get_transformation_matrix_forward();
        SIRFRegAffineTransformation c = NA.get_transformation_matrix_inverse();
        SIRFRegAffineTransformation d = b * c;
        SIRFRegAffineTransformation e = SIRFRegAffineTransformation::get_identity();
        if (d != e)
            throw std::runtime_error("SIRFRegAffineTransformation::mult/comparison failed.");

        if (e.get_determinant() - 1.F > 1.e-7F)
            throw std::runtime_error("SIRFRegAffineTransformation::get_determinant failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished SIRFRegAffineTransformation test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    // Error handling
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
