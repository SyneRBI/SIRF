/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

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
#include "sirf/Reg/NiftyAladinSym.h"
#include "sirf/Reg/NiftyF3dSym.h"
#include "sirf/Reg/NiftyResample.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/ImageWeightedMean.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/Quaternion.h"
#include <memory>
#include <numeric>
#ifdef SIRF_SPM
#include "sirf/Reg/SPMRegistration.h"
#endif

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
    const std::string save_nifti_image_upsample                  = output_prefix   + "save_NiftiImageData_upsample.nii";
    const std::string save_nifti_image_downsample                = output_prefix   + "save_NiftiImageData_downsample.nii";
    const std::string save_nifti_image_up_downsample             = output_prefix   + "save_NiftiImageData_upsample_downsample.nii";
    const std::string flo_aladin_as_unsigned_int                 = output_prefix   + "flo_aladin_as_unsigned_int.nii";
    const std::string aladin_warped            = output_prefix   + "aladin_warped.nii";
    const std::string f3d_warped               = output_prefix   + "f3d_warped.nii";
    const std::string TM_forward               = output_prefix   + "TM_forward.txt";
    const std::string TM_inverse               = output_prefix   + "TM_inverse.txt";
    const std::string aladin_def_forward       = output_prefix   + "aladin_def_forward.nii";
    const std::string aladin_def_inverse_xyz   = output_prefix   + "aladin_def_inverse_%s.nii";
    const std::string aladin_def_inverse       = output_prefix   + "aladin_def_inverse.nii";
    const std::string aladin_def_fwd_inv       = output_prefix   + "aladin_def_fwd_then_inv.nii";
    const std::string aladin_disp_forward      = output_prefix   + "aladin_disp_forward.nii";
    const std::string aladin_disp_inverse      = output_prefix   + "aladin_disp_inverse_%s.nii";
    const std::string f3d_disp_forward         = output_prefix   + "f3d_disp_forward.nii";
    const std::string f3d_disp_inverse         = output_prefix   + "f3d_disp_inverse_%s.nii";
    const std::string f3d_def_forward          = output_prefix   + "f3d_def_forward.nii";
    const std::string f3d_def_inverse          = output_prefix   + "f3d_def_inverse_%s.nii";
    const std::string rigid_resample           = output_prefix   + "rigid_resample.nii";
    const std::string nonrigid_resample_disp   = output_prefix   + "nonrigid_resample_disp.nii";
    const std::string nonrigid_resample_def    = output_prefix   + "nonrigid_resample_def.nii";
    const std::string niftymomo_resample_adj   = output_prefix   + "niftymomo_resample_adj.nii";
    const std::string output_weighted_mean     = output_prefix   + "weighted_mean.nii";
    const std::string output_weighted_mean_def = output_prefix   + "weighted_mean_def.nii";
    const std::string output_float             = output_prefix   + "reg_aladin_float.nii";
    const std::string spm_working_folder       = output_prefix   + "spm_working_folder";
    const std::string spm_working_folder2      = output_prefix   + "spm_working_folder2";

    const std::shared_ptr<const NiftiImageData3D<float> > ref_aladin(new NiftiImageData3D<float>( ref_aladin_filename ));
    const std::shared_ptr<const NiftiImageData3D<float> > flo_aladin(new NiftiImageData3D<float>( flo_aladin_filename ));
    const std::shared_ptr<const NiftiImageData3D<float> > ref_f3d   (new NiftiImageData3D<float>(   ref_f3d_filename  ));
    const std::shared_ptr<const NiftiImageData3D<float> > flo_f3d   (new NiftiImageData3D<float>(   flo_f3d_filename  ));

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // default constructor
        NiftiImageData<float> a;

        // Read from file
        NiftiImageData<float> b(ref_aladin_filename);

        // Save to file
        b.write(save_nifti_image);

        // Fill
        b.fill(100);

        // Get max
        if (fabs(b.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData fill()/get_max() failed.");

        // Get min
        if (fabs(b.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData fill()/get_min() failed.");

        // Deep copy
        NiftiImageData<float> d = b;
        std::cout << "\ndone.\n";
        if (d.get_raw_nifti_sptr() == b.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData deep_copy failed.");
        if (b != d)
            throw std::runtime_error("NiftiImageData deep_copy failed.");

        // Addition
        NiftiImageData3D<float> e = d + d;
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
        NiftiImageData3D<float> q = e + 1;
        if (fabs(q.get_max() - (e.get_max() + 1.F)) > 0.0001F)
            throw std::runtime_error("NiftiImageData __add__ val failed.");

        // Subtract num from image
        NiftiImageData3D<float> r = e - 1;
        if (fabs(r.get_max() - (e.get_max() - 1.F)) > 0.0001F)
            throw std::runtime_error("NiftiImageData __sub__ val failed.");

        // Multiply image by num
        NiftiImageData3D<float> s = e * 10;
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
        (*ref_aladin)(idx);

        // Test get_norm
        if (ref_aladin->get_norm(*flo_aladin) < 1.e-7F)
            throw std::runtime_error("NiftiImageData get_norm() failed.");

        // Test saving to datatype
        ref_aladin->write(output_float,NIFTI_TYPE_FLOAT32);
        NiftiImageData3D<float> ref_aladin_float(output_float);
        for (int i=0; i<int(ref_aladin->get_raw_nifti_sptr()->nvox); ++i)
            if (ref_aladin_float(i) - (*ref_aladin)(i) > 1.e-7F)
                throw std::runtime_error("NiftiImageData3D::write()/change_datatype() failed.");

        // Test print methods
        q.print_header();
        NiftiImageData<float>::print_headers({&q,&s});

        // Crop image
        int min[7], max[7];
        for (int i=0; i<7; ++i) {
            min[i] = 0;
            max[i] = f[i+1] - 1;
        }
        max[2] = 62;
        NiftiImageData<float> z = e;
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

        // Test creating an image from an array
        // Copy image, and convert float array to unsigned int array
        b = NiftiImageData<float>(flo_aladin_filename);
        unsigned int *data_array = new unsigned int[b.get_raw_nifti_sptr()->nvox];
        for (unsigned i=0; i<b.get_raw_nifti_sptr()->nvox; ++i)
            data_array[i] = static_cast<unsigned int>(b(i));
        // Construct image
        NiftiImageData<float> t(data_array, *b.get_geom_info_sptr());
        // Delete array
        delete [] data_array;
        data_array = nullptr;
        // Check values are still the same (there would be rounding involved, since we originally turn
        // the float into unsigned int and then back). However, this is ok, as we know that the input
        // image was already of type NIFTI_TYPE_UINT8 (unsigned char).
        if (b != t)
            throw std::runtime_error("NiftiImageData constructor from array.");
        // Save to file (useful for UI comparison)
        t.write(flo_aladin_as_unsigned_int);

        // Check upsampling/downsampling
        NiftiImageData<float> u(ref_aladin_filename);
        float *pixdim = u.get_raw_nifti_sptr()->pixdim;
        float original_spacing[3]    = {pixdim[1],       pixdim[2],       pixdim[3]};
        float upsampled_spacing[3]   = {pixdim[1] / 2.F, pixdim[2] / 4.F, pixdim[3]};
        float downsampled_spacing[3] = {pixdim[1] * 2.F, pixdim[2] * 4.F, pixdim[3]};
        // Downsample
        NiftiImageData<float> v = u;
        v.set_voxel_spacing(downsampled_spacing,3);
        v.write(save_nifti_image_downsample);
        // Upsample then downsample, check nothing has changed
        NiftiImageData<float> w = u;
        w.set_voxel_spacing(upsampled_spacing,0);
        w.write(save_nifti_image_upsample);
        NiftiImageData<float> x = w;
        x.set_voxel_spacing(original_spacing,0);
        x.write(save_nifti_image_up_downsample);
        NiftiImageData<float>::print_headers({&u, &v, &w, &x});
        if (x != u)
            throw std::runtime_error("NiftiImageData::upsample()/downsample() failed.");

        // Test inner product
        NiftiImageData<float> y = x;
        for (unsigned i=0; i<x.get_num_voxels(); ++i)
            x(int(i)) = static_cast<float>(i);
        for (unsigned i=0; i<x.get_num_voxels(); ++i)
            y(int(i)) = static_cast<float>(3*x.get_num_voxels()-i);
        const float inner = x.get_inner_product(y);

        // Do it with vectors to check
        const float *x_begin = &static_cast<const float*>(x.get_raw_nifti_sptr()->data)[0];
        const float *y_begin = &static_cast<const float*>(y.get_raw_nifti_sptr()->data)[0];
        const float *x_end   = &static_cast<const float*>(x.get_raw_nifti_sptr()->data)[0] + x.get_num_voxels();
        const float inner_vec = std::inner_product(x_begin, x_end, y_begin, 0.f);

        if (std::abs(inner-inner_vec) > 1e-4f)
            throw std::runtime_error("NiftiImageData::get_inner_product() failed.");

        // Test contains NaNs
        x.fill(0.f);
        if (x.get_contains_nans())
            throw std::runtime_error("NiftiImageData::get_contains_nans() 1 failed.");
        x(0) = NAN;
        if (!x.get_contains_nans())
            throw std::runtime_error("NiftiImageData::get_contains_nans() 2 failed.");

        // Test that eg im += 5 gives same as im = im + 5
        NiftiImageData<float> aa = *flo_aladin->clone();
        NiftiImageData<float> bb = *flo_aladin->clone();
        aa = aa + 5;
        bb += 5;
        if (bb != aa)
            throw std::runtime_error("NiftiImageData::+= (scalar) failed");
        aa = aa - 5;
        bb -= 5;
        if (bb != aa)
            throw std::runtime_error("NiftiImageData::-= (scalar) failed");
        aa = aa * 5;
        bb *= 5;
        if (bb != aa)
            throw std::runtime_error("NiftiImageData::*= failed");
        aa = aa / 5;
        bb /= 5;
        if (bb != aa)
            throw std::runtime_error("NiftiImageData::/= failed");

        aa = aa + aa;
        bb += bb;
        if (bb != aa)
            throw std::runtime_error("NiftiImageData::+= failed");
        aa = aa - aa;
        bb -= bb;
        if (bb != aa)
            throw std::runtime_error("NiftiImageData::-= failed");

        // Pad then crop, should be the same
        NiftiImageData<float> cc = aa;
        const int * const original_dims = aa.get_dimensions();
        int pad_in_min_dir[7] = { 1, 2, 3, 0, 0, 0, 0 };
        int pad_in_max_dir[7] = { 4, 5, 6, 0, 0, 0, 0 };
        cc.pad(pad_in_min_dir,pad_in_max_dir, 100.f);
        const int * const padded_dims = cc.get_dimensions();
        for (unsigned i=0; i<7; ++i)
            if (padded_dims[i+1] != original_dims[i+1] + pad_in_min_dir[i] + pad_in_max_dir[i])
                throw std::runtime_error("NiftiImageData::pad failed");
        // Crop back to beginning
        int cropped_min_dir[7], cropped_max_dir[7];
        for (unsigned i=0; i<7; ++i) {
            cropped_min_dir[i] = pad_in_min_dir[i];
            cropped_max_dir[i] = original_dims[i+1] + cropped_min_dir[i] - 1;
        }
        cc.crop(cropped_min_dir, cropped_max_dir);
        if (aa != cc)
            throw std::runtime_error("NiftiImageData::pad/crop failed");


        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData3D test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // default constructor
        NiftiImageData3D<float> a;

        // Read from file
        NiftiImageData3D<float> b(ref_aladin_filename);

        // Save to file
        b.write(save_nifti_image_3d);

        // Fill
        b.fill(100);

        // Get max
        if (fabs(b.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3D fill()/get_max() failed.");

        // Get min
        if (fabs(b.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3D fill()/get_min() failed.");

        // Deep copy
        NiftiImageData3D<float> d = b;
        if (d.get_raw_nifti_sptr() == b.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData3D deep_copy failed.");
        if (d != b)
            throw std::runtime_error("NiftiImageData3D deep_copy failed.");

        // Addition
        NiftiImageData3D<float> e = d + d;
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

        // Test creating an image from an array
        // Copy image, and convert float array to unsigned int array
        b = NiftiImageData3D<float>(flo_aladin_filename);
        unsigned int *data_array = new unsigned int[b.get_raw_nifti_sptr()->nvox];
        for (unsigned i=0; i<b.get_raw_nifti_sptr()->nvox; ++i)
            data_array[i] = static_cast<unsigned int>(b(i));
        // Construct image
        NiftiImageData3D<float> t(data_array, *b.get_geom_info_sptr());
        // Delete array
        delete [] data_array;
        data_array = nullptr;
        // Check values are still the same (there would be rounding involved, since we originally turn
        // the float into unsigned int and then back). However, this is ok, as we know that the input
        // image was already of type NIFTI_TYPE_UINT8 (unsigned char).
        if (b != t)
            throw std::runtime_error("NiftiImageData3D constructor from array.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3D test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData3DTensor test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImageData3DTensor from NiftiImageData3D
        NiftiImageData3DTensor<float> b;
        b.create_from_3D_image(*ref_aladin);

        // Save to file
        b.write(save_nifti_image_3d_tensor_not_split);
        b.write_split_xyz_components(save_nifti_image_3d_tensor_split);

        // Constructor from file
        NiftiImageData3DTensor<float> c(save_nifti_image_3d_tensor_not_split);

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DTensor fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DTensor fill()/get_min() failed.");

        // Deep copy
        NiftiImageData3DTensor<float> d = c;
        if (d.get_raw_nifti_sptr() == c.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData3DTensor deep_copy failed.");
        if (d != c)
            throw std::runtime_error("NiftiImageData3DTensor deep_copy failed.");

        // Addition
        NiftiImageData3DTensor<float> e = d + d;
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
        NiftiImageData3D<float> im1 = *ref_aladin;
        NiftiImageData3D<float> im2 = *ref_aladin;
        NiftiImageData3D<float> im3 = *ref_aladin;
        im1.fill(30.F);
        im2.fill(20.F);
        im3.fill(-10.F);
        std::vector<NiftiImageData3D<float> > ims = {im1, im2, im3};
        NiftiImageData3DTensor<float> h(im1, im2, im3);
        for (int i=0; i<3; ++i)
            if (*h.get_tensor_component(i) != ims.at(i))
                throw std::runtime_error("NiftiImageData3DTensor 3ims->tensor->3ims failed on idx " + std::to_string(i) + ".");

        // Test flip components
        h.flip_component(0);
        if (fabs(h.get_max() - 20.F) > 1.e-7F )
            throw std::runtime_error("NiftiImageData3DTensor flip_component() failed.");
        if (fabs(h.get_min() + 30.F) > 1.e-7F )
            throw std::runtime_error("NiftiImageData3DTensor flip_component() failed.");

        // Test creating an image from an array
        float *data_array = new float[h.get_raw_nifti_sptr()->nvox];
        for (unsigned i=0; i<h.get_raw_nifti_sptr()->nvox; ++i)
            data_array[i] = static_cast<float>(h(i));
        // Construct image
        NiftiImageData3DTensor<float> t(data_array, *h.get_geom_info_sptr());
        // Delete array
        delete [] data_array;
        data_array = nullptr;
        if (h != t)
            throw std::runtime_error("NiftiImageData3DTensor constructor from array.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3DTensor test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData3DDisplacement test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImageData3DDisplacement from NiftiImageData3D
        NiftiImageData3DDisplacement<float> b;
        b.create_from_3D_image(*ref_aladin);

        // Save to file
        b.write(save_nifti_image_3d_displacement_not_split);
        b.write_split_xyz_components(save_nifti_image_3d_displacement_split);

        // Constructor from file
        NiftiImageData3DDisplacement<float> c(save_nifti_image_3d_displacement_not_split);

        // Constructor from tensor
        NiftiImageData3DTensor<float> x(save_nifti_image_3d_displacement_not_split);
        NiftiImageData3DDisplacement<float> y(x);

        // Constructor from general
        NiftiImageData<float> q(save_nifti_image_3d_displacement_not_split);
        NiftiImageData3DDisplacement<float> r(q);

        // Constructor from single components
        NiftiImageData3DDisplacement<float> h(*ref_aladin,*ref_aladin,*ref_aladin);
        for (int i=0; i<3; ++i)
            if (*h.get_tensor_component(i) != *ref_aladin)
                throw std::runtime_error("NiftiImageData3DDisplacement 3ims->tensor->3ims failed on idx " + std::to_string(i) + ".");

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DDisplacement fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DDisplacement fill()/get_min() failed.");

        // Deep copy
        NiftiImageData3DDisplacement<float> d = c;
        if (d.get_raw_nifti_sptr() == c.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData3DDisplacement deep_copy failed.");
        if (d != c)
            throw std::runtime_error("NiftiImageData3DDisplacement deep_copy failed.");

        // Addition
        NiftiImageData3DDisplacement<float> e = d + d;
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

        // Test creating an image from an array
        float *data_array = new float[h.get_raw_nifti_sptr()->nvox];
        for (unsigned i=0; i<h.get_raw_nifti_sptr()->nvox; ++i)
            data_array[i] = static_cast<float>(h(i));
        // Construct image
        NiftiImageData3DDisplacement<float> t(data_array, *h.get_geom_info_sptr());
        // Delete array
        delete [] data_array;
        data_array = nullptr;
        if (h != t)
            throw std::runtime_error("NiftiImageData3DDisplacement constructor from array.");

        // Check upsampling/downsampling
        NiftiImageData3DDisplacement<float> u(save_nifti_image_3d_displacement_not_split);
        float *pixdim = u.get_raw_nifti_sptr()->pixdim;
        float original_spacing[3]    = {pixdim[1],       pixdim[2],       pixdim[3]};
        float upsampled_spacing[3]   = {pixdim[1] / 2.F, pixdim[2] / 4.F, pixdim[3]};
        float downsampled_spacing[3] = {pixdim[1] * 2.F, pixdim[2] * 4.F, pixdim[3]};
        // Downsample
        NiftiImageData3DDisplacement<float> v = u;
        v.set_voxel_spacing(downsampled_spacing,3);
        // Upsample then downsample, check nothing has changed
        NiftiImageData3DDisplacement<float> w = u;
        w.set_voxel_spacing(upsampled_spacing,0);
        NiftiImageData3DDisplacement<float> z = w;
        z.set_voxel_spacing(original_spacing,0);
        NiftiImageData3DDisplacement<float>::print_headers({&u, &v, &w, &z});
        if (z != u)
            throw std::runtime_error("NiftiImageData3DDisplacement::upsample()/downsample() failed.");


        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3DDisplacement test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftiImageData3DDeformation test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Create NiftiImageData3DDeformation from NiftiImageData3D
        NiftiImageData3DDeformation<float> b;
        b.create_from_3D_image(*ref_aladin);

        // Save to file
        b.write(save_nifti_image_3d_deformation_not_split);
        b.write_split_xyz_components(save_nifti_image_3d_deformation_split);

        // Constructor from file
        NiftiImageData3DDeformation<float> c(save_nifti_image_3d_deformation_not_split);

        // Constructor from tensor
        NiftiImageData3DTensor<float> x(save_nifti_image_3d_deformation_not_split);
        NiftiImageData3DDeformation<float> y(x);

        // Constructor from general
        NiftiImageData<float> q(save_nifti_image_3d_deformation_not_split);
        NiftiImageData3DDeformation<float> r(q);

        // Constructor from single components
        NiftiImageData3DDeformation<float> h(*ref_aladin,*ref_aladin,*ref_aladin);
        for (int i=0; i<3; ++i)
            if (*h.get_tensor_component(i) != *ref_aladin)
                throw std::runtime_error("NiftiImageData3DDeformation 3ims->tensor->3ims failed on idx " + std::to_string(i) + ".");

        // Fill
        c.fill(100);

        // Get max
        if (fabs(c.get_max() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DDeformation fill()/get_max() failed.");

        // Get min
        if (fabs(c.get_min() - 100) > 1.e-5F)
            throw std::runtime_error("NiftiImageData3DDeformation fill()/get_min() failed.");

        // Deep copy
        NiftiImageData3DDeformation<float> d = c;
        if (d.get_raw_nifti_sptr() == c.get_raw_nifti_sptr())
            throw std::runtime_error("NiftiImageData3DDeformation deep_copy failed.");
        if (d != c)
            throw std::runtime_error("NiftiImageData3DDeformation deep_copy failed.");

        // Addition
        NiftiImageData3DDeformation<float> e = d + d;
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

        // Test creating an image from an array
        float *data_array = new float[h.get_raw_nifti_sptr()->nvox];
        for (unsigned i=0; i<h.get_raw_nifti_sptr()->nvox; ++i)
            data_array[i] = static_cast<float>(h(i));
        // Construct image
        NiftiImageData3DDeformation<float> t(data_array, *h.get_geom_info_sptr());
        // Delete array
        delete [] data_array;
        data_array = nullptr;
        if (h != t)
            throw std::runtime_error("NiftiImageData3DDeformation constructor from array.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3DDeformation test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    NiftyAladinSym<float> NA;
    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Nifty aladin test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Print all available methods
        NiftyAladinSym<float>::print_all_wrapped_methods();

        // First set up some masks
        std::shared_ptr<NiftiImageData3D<float> > ref_mask = ref_aladin->clone();
        std::shared_ptr<NiftiImageData3D<float> > flo_mask = flo_aladin->clone();
        ref_mask->fill(1.F);
        flo_mask->fill(1.F);

        NA.set_reference_image               (            ref_aladin         );
        NA.set_floating_image                (            flo_aladin         );
        NA.set_parameter_file                (      parameter_file_aladin    );
        NA.set_parameter("SetInterpolationToCubic");
        NA.set_parameter("SetLevelsToPerform","1");
        NA.set_parameter("SetMaxIterations","5");
        NA.set_parameter("SetPerformRigid","1");
        NA.set_parameter("SetPerformAffine","0");
        NA.set_reference_mask(ref_mask);
        NA.set_floating_mask(flo_mask);
        NA.process();

        // Get outputs
        const NiftiImageData3D<float>               warped_     = dynamic_cast<const NiftiImageData3D<float>&>(*NA.get_output_sptr());
        const AffineTransformation<float>         TM_forward_   = dynamic_cast<const AffineTransformation<float>&> (*NA.get_transformation_matrix_forward_sptr());
        const AffineTransformation<float>         TM_inverse_   = dynamic_cast<const AffineTransformation<float>&> (*NA.get_transformation_matrix_forward_sptr());
        const NiftiImageData3DDeformation<float>  def_forward_  = dynamic_cast<const NiftiImageData3DDeformation<float>&> (*NA.get_deformation_field_forward_sptr());
        const NiftiImageData3DDeformation<float>  def_inverse_  = dynamic_cast<const NiftiImageData3DDeformation<float>&> (*NA.get_deformation_field_inverse_sptr());
        const NiftiImageData3DDisplacement<float> disp_forward_ = dynamic_cast<const NiftiImageData3DDisplacement<float>&>(*NA.get_displacement_field_forward_sptr());
        const NiftiImageData3DDisplacement<float> disp_inverse_ = dynamic_cast<const NiftiImageData3DDisplacement<float>&>(*NA.get_displacement_field_inverse_sptr());

        // Check registration with filenames
        NA.set_reference_image_filename(ref_aladin_filename);
        NA.set_floating_image_filename(flo_aladin_filename);
        NA.process();

        const std::shared_ptr<const NiftiImageData3D<float> >             warped_sptr       = std::dynamic_pointer_cast<const NiftiImageData3D<float> >(NA.get_output_sptr());
        const std::shared_ptr<const AffineTransformation<float> >         TM_forward_sptr   = std::dynamic_pointer_cast<const AffineTransformation<float> > (NA.get_transformation_matrix_forward_sptr());
        const std::shared_ptr<const AffineTransformation<float> >         TM_inverse_sptr   = std::dynamic_pointer_cast<const AffineTransformation<float> > (NA.get_transformation_matrix_forward_sptr());
        const std::shared_ptr<const NiftiImageData3DDeformation<float> >  def_forward_sptr  = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> > (NA.get_deformation_field_forward_sptr());
        const std::shared_ptr<const NiftiImageData3DDeformation<float> >  def_inverse_sptr  = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> > (NA.get_deformation_field_inverse_sptr());
        const std::shared_ptr<const NiftiImageData3DDisplacement<float> > disp_forward_sptr = std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(NA.get_displacement_field_forward_sptr());
        const std::shared_ptr<const NiftiImageData3DDisplacement<float> > disp_inverse_sptr = std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(NA.get_displacement_field_inverse_sptr());

        if (*warped_sptr           != warped_       ||
                *TM_forward_sptr   != TM_forward_   ||
                *TM_inverse_sptr   != TM_inverse_   ||
                *def_forward_sptr  != def_forward_  ||
                *def_inverse_sptr  != def_inverse_  ||
                *disp_inverse_sptr != disp_inverse_ ||
                *disp_inverse_sptr != disp_inverse_)
            throw std::runtime_error("Error doing registration via filename");

        warped_sptr->write    (      aladin_warped    );
        TM_forward_sptr->write(       TM_forward      );
        TM_inverse_sptr->write(       TM_inverse      );
        disp_forward_sptr->write(aladin_disp_forward);
        disp_inverse_sptr->write_split_xyz_components(aladin_disp_inverse);
        def_forward_sptr->write(aladin_def_forward);
        def_inverse_sptr->write_split_xyz_components(aladin_def_inverse_xyz);
        def_inverse_sptr->write(aladin_def_inverse);

        // forward TM
        TM_forward_sptr->print();

        // Inverse TM
        const AffineTransformation<float> &inverse_tm = *NA.get_transformation_matrix_inverse_sptr();
        inverse_tm.print();

        // Test converting disp to def
        NiftiImageData3DDeformation<float> a(*disp_forward_sptr);
        if (a != *def_forward_sptr)
            throw std::runtime_error("NiftiImageData3DDeformation::create_from_disp() failed.");

        // Test converting def to disp
        NiftiImageData3DDisplacement<float> b(*def_forward_sptr);
        if (b != *disp_forward_sptr)
            throw std::runtime_error("NiftiImageData3DDisplacement::create_from_def() failed.");

        // Check NiftiImageData3DDeformation::get_inverse()
        const std::shared_ptr<const NiftiImageData3DDeformation<float> > def_fwd_then_inv_sptr =
                def_forward_sptr->get_inverse(flo_aladin);
        def_fwd_then_inv_sptr->write(aladin_def_fwd_inv);
        NiftiImageData<float>::print_headers({&*ref_aladin, &*flo_aladin, &*def_inverse_sptr, &*def_fwd_then_inv_sptr});

        // Reference forward with def_inv
        NiftyResample<float> resample;
        resample.set_reference_image(flo_aladin);
        resample.set_floating_image(ref_aladin);
        resample.set_padding_value(0.f);
        resample.set_interpolation_type_to_linear();
        resample.add_transformation(def_inverse_sptr);
        const std::shared_ptr<const NiftiImageData<float> > out1_sptr = std::dynamic_pointer_cast<const NiftiImageData<float> >(resample.forward(ref_aladin));

        // Reference forward with def_fwd_then_inv_sptr
        resample.clear_transformations();
        resample.add_transformation(def_fwd_then_inv_sptr);
        const std::shared_ptr<const NiftiImageData<float> > out2_sptr = std::dynamic_pointer_cast<const NiftiImageData<float> >(resample.forward(ref_aladin));

        NiftiImageData<float>::print_headers({&*out1_sptr, &*out2_sptr});

        if (*out1_sptr != *out2_sptr)
            throw std::runtime_error("NiftiImageData3DDeformation::get_inverse() failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished Nifty aladin test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Nifty f3d test..\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Print all available methods
        NiftyF3dSym<float>::print_all_wrapped_methods();

        // First set up some masks
        std::shared_ptr<NiftiImageData3D<float> > ref_mask = ref_f3d->clone();
        std::shared_ptr<NiftiImageData3D<float> > flo_mask = flo_f3d->clone();
        ref_mask->fill(1.F);
        flo_mask->fill(1.F);

        NiftyF3dSym<float> NF;
        NF.set_reference_image               (           ref_f3d          );
        NF.set_floating_image                (           flo_f3d          );
        NF.set_parameter_file                (     parameter_file_f3d     );
        NF.set_reference_time_point          (             1              );
        NF.set_floating_time_point           (             1              );
        NF.set_reference_mask(ref_mask);
        NF.set_floating_mask(flo_mask);
        NF.process();

        // Get outputs
        std::shared_ptr<const NiftiImageData3D<float> >             warped_sptr       = std::dynamic_pointer_cast<const NiftiImageData3D<float> >            (NF.get_output_sptr());
        std::shared_ptr<const NiftiImageData3DDeformation<float> >  def_forward_sptr  = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> > (NF.get_deformation_field_forward_sptr());
        std::shared_ptr<const NiftiImageData3DDeformation<float> >  def_inverse_sptr  = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> > (NF.get_deformation_field_inverse_sptr());
        std::shared_ptr<const NiftiImageData3DDisplacement<float> > disp_forward_sptr = std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(NF.get_displacement_field_forward_sptr());
        std::shared_ptr<const NiftiImageData3DDisplacement<float> > disp_inverse_sptr = std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(NF.get_displacement_field_inverse_sptr());

        warped_sptr->write      (  f3d_warped   );
        def_forward_sptr->write (f3d_def_forward);
        def_inverse_sptr->write_split_xyz_components(f3d_def_inverse);
        disp_forward_sptr->write(f3d_disp_forward);
        disp_inverse_sptr->write_split_xyz_components(f3d_disp_inverse);


        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished Nifty f3d test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting transformations test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Get outputs from aladin
        std::shared_ptr<const NiftiImageData3DDeformation<float> >  def_forward_sptr  = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> > (NA.get_deformation_field_forward_sptr());
        std::shared_ptr<const NiftiImageData3DDisplacement<float> > disp_forward_sptr = std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(NA.get_displacement_field_forward_sptr());

        // Affine
        std::cout << "\nTesting affine...\n";
        AffineTransformation<float> a1;
        AffineTransformation<float> a2(TM_forward);
        AffineTransformation<float> a3(*NA.get_transformation_matrix_forward_sptr());

        // Displacement
        std::cout << "\nTesting displacement...\n";
        const NiftiImageData3DDisplacement<float> b3(*disp_forward_sptr);

        // Deformation
        std::cout << "\nTesting deformation...\n";
        NiftiImageData3DDeformation<float> c3(*def_forward_sptr);

        // Get as deformations
        NiftiImageData3DDeformation<float> a_def = a3.get_as_deformation_field(*ref_aladin);
        NiftiImageData3DDeformation<float> b_def = b3.get_as_deformation_field(*ref_aladin);
        NiftiImageData3DDeformation<float> c_def = c3.get_as_deformation_field(*ref_aladin);
        if (a_def != *def_forward_sptr)
            throw std::runtime_error("AffineTransformation::get_as_deformation_field failed.");
        if (b_def != *def_forward_sptr)
            throw std::runtime_error("NiftiImageData3DDisplacement::get_as_deformation_field failed.");
        if (c_def != *def_forward_sptr)
            throw std::runtime_error("NiftiImageData3DDeformation::get_as_deformation_field failed.");

        // Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
        AffineTransformation<float> tm_iden;
        AffineTransformation<float> trans_aff_iden(tm_iden);
        std::vector<std::shared_ptr<const Transformation<float> > > vec;
        vec.push_back(std::make_shared<const AffineTransformation<float> >(trans_aff_iden));
        vec.push_back(std::make_shared<const AffineTransformation<float> >(trans_aff_iden));
        vec.push_back(std::make_shared<const NiftiImageData3DDeformation<float> >(c3));
        NiftiImageData3DDeformation<float> composed =
                NiftiImageData3DDeformation<float>::compose_single_deformation(vec, *ref_aladin);
        if (composed.get_as_deformation_field(*ref_aladin) != *def_forward_sptr)
            throw std::runtime_error("NiftiImageData3DDeformation::compose_single_deformation failed.");

        // Test get_inverse
        AffineTransformation<float> tm_inv = tm_iden.get_inverse();

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished transformations test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Nifty resample test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        std::shared_ptr<const Transformation<float> > tm_iden  = std::make_shared<const AffineTransformation<float> >();
        std::shared_ptr<const Transformation<float> > tm       = NA.get_transformation_matrix_forward_sptr();
        std::shared_ptr<const Transformation<float> > disp     = NA.get_displacement_field_forward_sptr();
        std::shared_ptr<const Transformation<float> > deff     = NA.get_deformation_field_forward_sptr();
        float padding_value = -20.f;

        std::cout << "Testing rigid resample...\n";
        NiftyResample<float> nr1;
        nr1.set_reference_image(ref_aladin);
        nr1.set_floating_image(flo_aladin);
        nr1.set_interpolation_type_to_cubic_spline(); // try different interpolations
        nr1.set_interpolation_type(NiftyResample<float>::CUBICSPLINE); // try different interpolations (cubic)
        nr1.add_transformation(tm_iden);
        nr1.clear_transformations();
        nr1.add_transformation(tm_iden);
        nr1.add_transformation(tm);
        nr1.process();
        nr1.get_output_sptr()->write(rigid_resample);

        std::cout << "Testing non-rigid displacement...\n";
        NiftyResample<float> nr2;
        nr2.set_reference_image(ref_aladin);
        nr2.set_floating_image(flo_aladin);
        nr2.set_interpolation_type_to_sinc(); // try different interpolations
        nr2.set_interpolation_type_to_linear(); // try different interpolations
        nr2.add_transformation(disp);
        nr2.set_padding_value(padding_value);
        nr2.process();
        const std::shared_ptr<const NiftiImageData<float> > nr2_output =
                std::dynamic_pointer_cast<const NiftiImageData<float> >(
                    nr2.get_output_sptr());
        nr2_output->write(nonrigid_resample_disp);

        if (std::abs(nr2_output->get_min() - padding_value) > 1e-4f) // only get exact value with linear inerpolation
            throw std::runtime_error("NiftyResample::set_padding_value failed.");

        std::cout << "Testing non-rigid deformation...\n";
        NiftyResample<float> nr3;
        nr3.set_reference_image(ref_aladin);
        nr3.set_floating_image(flo_aladin);
        nr3.set_interpolation_type_to_nearest_neighbour(); // try different interpolations
        nr3.add_transformation(deff);
        nr3.set_interpolation_type_to_linear();
        nr3.process();
        nr3.get_output_sptr()->write(nonrigid_resample_def);

        // Check that the following give the same result
        //      out = resample.forward(in)
        //      resample.forward(out, in)
        const std::shared_ptr<const NiftiImageData<float> > out1_sptr =
                std::dynamic_pointer_cast<const NiftiImageData<float> >(
                    nr3.forward(flo_aladin));

        const std::shared_ptr<NiftiImageData<float> > out2_sptr = ref_aladin->clone();
        nr3.forward(out2_sptr, flo_aladin);

        if (*out1_sptr != *out2_sptr)
            throw std::runtime_error("out = NiftyResample::forward(in) and NiftyResample::forward(out, in) do not give same result.");

        // TODO this doesn't work. For some reason (even with NiftyReg directly), resampling with the TM from the registration
        // doesn't give the same result as the output from the registration itself (even with same interpolations). Even though
        // ref and flo images are positive, the output of the registration can be negative. This implies that linear interpolation
        // is not used to generate final image. You would hope it's used throughout the registration process, otherwise why is it there?
        // i.e., reg_aladin != reg_resample when reg_resample uses the transformation matrix from reg_aladin
        /*if (NA.get_output() != nr1.get_output())
            throw std::runtime_error("compose_transformations_into_single_deformation failed.");*/

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished Nifty resample test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting NiftyMoMo test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // The forward and the adjoint should meet the following criterion:
        //      |<x, Ty> - <y, Tsx>| / 0.5*(|<x, Ty>|+|<y, Tsx>|) < epsilon
        // for all images x and y, where T is the transform and Ts is the adjoint.

        const std::shared_ptr<const NiftiImageData<float> > x =
                std::make_shared<const NiftiImageData<float> >(ref_aladin_filename);
        const std::shared_ptr<AffineTransformation<float> > T =
                std::make_shared<AffineTransformation<float> >(*
                NA.get_transformation_matrix_forward_sptr());
        const std::shared_ptr<NiftiImageData<float> > y  =
                std::make_shared<NiftiImageData3D<float> >(flo_aladin_filename);

        // Add in a magnification to make things interesting
        (*T)[0][0] = 1.5f;

        // make it slightly unsquare to spice things up
        int min_idx[7] = {0,1,2,-1,-1,-1,-1};
        const int *y_dims = y->get_dimensions();
        int max_idx[7] = {y_dims[1]-3,y_dims[2]-1,y_dims[3]-5-1,-1,-1,-1};
        y->crop(min_idx,max_idx);

        NiftyResample<float> nr;
        nr.set_reference_image(x);
        nr.set_floating_image(y);
        nr.set_interpolation_type(Resample<float>::LINEAR);
        nr.add_transformation(T);

        // Do the forward
        const std::shared_ptr<const NiftiImageData<float> > Ty =
                std::dynamic_pointer_cast<const NiftiImageData<float> >(
                    nr.forward(y));

        // Do the adjoint
        const std::shared_ptr<const NiftiImageData<float> > Tsx =
                std::dynamic_pointer_cast<const NiftiImageData<float> >(
                    nr.adjoint(x));

        // Check the adjoint is truly the adjoint with: |<x, Ty> - <y, Tsx>| / 0.5*(|<x, Ty>|+|<y, Tsx>|) < epsilon
        float inner_x_Ty  = x->get_inner_product(*Ty);
        float inner_y_Tsx = y->get_inner_product(*Tsx);
        float adjoint_test = std::abs(inner_x_Ty - inner_y_Tsx) / (0.5f * (std::abs(inner_x_Ty) +std::abs(inner_y_Tsx)));
        std::cout << "\n<x, Ty>  = " << inner_x_Ty << "\n";
        std::cout << "<y, Tsx> = " << inner_y_Tsx << "\n";
        std::cout << "|<x, Ty> - <y, Tsx>| / 0.5*(|<x, Ty>|+|<y, Tsx>|) = " << adjoint_test << "\n";
        if (adjoint_test > 1e-4F)
            throw std::runtime_error("NiftyResample::adjoint() failed");

        // Check that the following give the same result
        //      out = resample.adjoint(in)
        //      resample.adjoint(out, in)
        const std::shared_ptr<const NiftiImageData<float> > out1_sptr =
                std::dynamic_pointer_cast<const NiftiImageData<float> >(
                    nr.adjoint(x));

        const std::shared_ptr<NiftiImageData<float> > out2_sptr = y->clone();
        nr.backward(out2_sptr, x);

        if (*out1_sptr != *out2_sptr)
            throw std::runtime_error("out = NiftyResample::adjoint(in) and NiftyResample::adjoint(out, in) do not give same result.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftyMoMo test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting weighted mean test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        //  Do 3D
        ImageWeightedMean<float> wm1;
        NiftiImageData3D<float> im1 = *ref_aladin;
        NiftiImageData3D<float> im2 = *ref_aladin;
        NiftiImageData3D<float> im3 = *ref_aladin;
        NiftiImageData3D<float> im4 = *ref_aladin;
        im1.fill(1);
        im2.fill(4);
        im3.fill(7);
        im4.fill(6);
        wm1.add_image(im1, 2.F);
        wm1.add_image(im2, 4.F);
        wm1.add_image(im3, 3.F);
        wm1.add_image(im4, 1.F);
        wm1.process();
        wm1.get_output_sptr()->write(output_weighted_mean);
        //  Answer should be 4.5, so compare it to that!
        NiftiImageData3D<float> res = *ref_aladin;
        res.fill(4.5F);

        if (*wm1.get_output_sptr() != res)
            throw std::runtime_error("ImageWeightedMean3D failed.");

        // Get def from aladin
        std::shared_ptr<const NiftiImageData3DDeformation<float> >  def_forward_sptr  = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> > (NA.get_deformation_field_forward_sptr());

        //  Do 4D
        ImageWeightedMean<float> wm2;
        NiftiImageData3DTensor<float> im4D1 = *def_forward_sptr;
        NiftiImageData3DTensor<float> im4D2 = *def_forward_sptr;
        NiftiImageData3DTensor<float> im4D3 = *def_forward_sptr;
        NiftiImageData3DTensor<float> im4D4 = *def_forward_sptr;
        im4D1.fill(1.F);
        im4D2.fill(4.F);
        im4D3.fill(7.F);
        im4D4.fill(6.F);
        wm2.add_image(im4D1, 2.F);
        wm2.add_image(im4D2, 4.F);
        wm2.add_image(im4D3, 3.F);
        wm2.add_image(im4D4, 1.F);
        wm2.process();
        wm2.get_output_sptr()->write(output_weighted_mean_def);
        //  Answer should be 4.5, so compare it to that!
        NiftiImageData3DTensor<float> res4D = *def_forward_sptr;
        res4D.fill(4.5);

        if (*wm2.get_output_sptr() != res4D)
            throw std::runtime_error("ImageWeightedMean3DTensor failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished weighted mean test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }
/* TODO UNCOMMENT WHEN GEOMETRICAL INFO IS IMPLEMENTED
    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting STIR to Nifti test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

            // Open stir image
            sirf::PETImageData pet_image_data(ref_aladin_filename);
            NiftiImageData3D<float> image_data_from_stir(pet_image_data);

            // Now fill the stir and nifti images with 1 and 100, respectively
            pet_image_data.fill(1.F);
            image_data_from_stir.fill(100.F);

            if (fabs(pet_image_data.data_sptr()->find_max() - image_data_from_stir.get_max()) < 1.e-5F)
                throw std::runtime_error("STIR & Nifti seem to share the same data pointers (their values should be different, but they're the same).");

            // Fill the stir image with the Nifti
            image_data_from_stir.copy_data_to(pet_image_data);
            if (fabs(pet_image_data.data_sptr()->find_max() - image_data_from_stir.get_max()) > 1.e-5F)
                throw std::runtime_error("NiftiImageData3D::copy_data_to failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished STIR to Nifti test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }
*/
    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting AffineTransformation test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Construct from file
        AffineTransformation<float> a(TM_forward);

        // Multiply forward and inverse, should equal identity
        AffineTransformation<float> b = *NA.get_transformation_matrix_forward_sptr();
        AffineTransformation<float> c = *NA.get_transformation_matrix_inverse_sptr();
        AffineTransformation<float> d = b * c;
        AffineTransformation<float> e;
        if (d != e)
            throw std::runtime_error("AffineTransformation::mult/comparison failed.");

        if (e.get_determinant() - 1.F > 1.e-7F)
            throw std::runtime_error("AffineTransformation::get_determinant failed.");

        // Test get_Euler_angles
        AffineTransformation<float> test_Eul;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                test_Eul[i][j]=0.F;
        test_Eul[0][2] =  1.F;
        test_Eul[1][1] = -1.F;
        test_Eul[2][0] = -1.F;
        test_Eul[3][3] =  1.F;
        // Example given by rotm2eul for MATLAB is [0 0 1; 0 -1 0; -1 0 0] -> XYZ = [-3.1416 1.5708 0]
        std::array<float,3> Eul = test_Eul.get_Euler_angles();
        std::array<float,3> Eul_expected{-3.1416F, 1.5708F, 0.F};
        for (unsigned i=0; i<3; ++i)
            if (std::abs(Eul[i] - Eul_expected[i]) > 1e-4F)
                throw std::runtime_error("AffineTransformation::get_Euler_angles failed.");

        // Average!
        Quaternion<float> quat_1(0.92707F,  0.02149F,   0.19191F,  0.32132F);
        Quaternion<float> quat_2(0.90361F,  0.0025836F, 0.097279F, 0.41716F);
        Quaternion<float> quat_3(0.75868F, -0.21289F,   0.53263F,  0.30884F);
        AffineTransformation<float> tm_1({0.F,0.F,0.F},quat_1);
        AffineTransformation<float> tm_2({0.F,0.F,0.F},quat_2);
        AffineTransformation<float> tm_3({0.F,0.F,0.F},quat_3);
        AffineTransformation<float> average = AffineTransformation<float>::get_average({tm_1,tm_2,tm_3});
        AffineTransformation<float> exptd_average;
        exptd_average[0][0] =  0.5836F;
        exptd_average[0][1] = -0.6736F;
        exptd_average[0][2] =  0.4535F;
        exptd_average[1][0] =  0.6007F;
        exptd_average[1][1] =  0.7339F;
        exptd_average[1][2] =  0.3171F;
        exptd_average[2][0] = -0.5464F;
        exptd_average[2][1] =  0.0874F;
        exptd_average[2][2] =  0.8329F;
        if (average != exptd_average)
            throw std::runtime_error("AffineTransformation::get_average() failed.");


        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished AffineTransformation test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Quaternion test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        // Construct TM
        AffineTransformation<float> rotm;
        for (unsigned i=0; i<4; ++i)
            for (unsigned j=0; j<4; ++j)
                rotm[i][j]=0.F;
        rotm[0][2] =  1.f;
        rotm[1][1] =  1.f;
        rotm[2][0] = -1.f;
        rotm[3][3] =  1.f;

        // Convert to quaternion
        Quaternion<float> quat(rotm);
        // Compare to expected values
        Quaternion<float> expt(0.707107f, 0.f, 0.707107f, 0.f);

        if (quat != expt)
            throw std::runtime_error("Quaternion from TM failed.");

        // Convert back to TM
        std::array<float,3> trans{0.F,0.F,0.F};
        AffineTransformation<float> affine(trans,quat);
        if (affine != rotm)
            throw std::runtime_error("TM to quaternion failed.");

        // Convert TM to quaternion
        Quaternion<float> quat2 = affine.get_quaternion();
        if (std::abs(quat.dot(quat2)) -1.f > 1.e-4f)
            throw std::runtime_error("AffineTransformation.get_quaternion()/Quaternion::dot() failed.");

        // Average!
        Quaternion<float> quat_1(0.92707F,  0.02149F,   0.19191F,  0.32132F);
        Quaternion<float> quat_2(0.90361F,  0.0025836F, 0.097279F, 0.41716F);
        Quaternion<float> quat_3(0.75868F, -0.21289F,   0.53263F,  0.30884F);
        Quaternion<float> average = Quaternion<float>::get_average({quat_1,quat_2,quat_3});
        Quaternion<float> exptd_average(0.88748F, -0.0647152F, 0.281671F, 0.35896F);

        if (average != exptd_average)
            throw std::runtime_error("Quaternion::get_average() failed.");

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished Quaternion test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }
#ifdef SIRF_SPM
    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting SPM test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        {

            // Resample an image with NiftyResample. Register SPM, check the result

            // TM
            std::array<float,3> translations = {5.f, 4.f, -5.f};
            std::array<float,3> euler_angles = {5.f, -2.f,  3.f};
            const std::shared_ptr<const AffineTransformation<float> > tm_sptr =
                    std::make_shared<const AffineTransformation<float> >(translations,euler_angles,true);

            NiftyResample<float> niftyreg_resampler;
            niftyreg_resampler.set_padding_value(0.f);
            niftyreg_resampler.set_reference_image(ref_aladin);
            niftyreg_resampler.set_floating_image(ref_aladin);
            niftyreg_resampler.add_transformation(tm_sptr);
            niftyreg_resampler.set_interpolation_type_to_linear();
            const std::shared_ptr<const ImageData> floating_sptr = niftyreg_resampler.forward(ref_aladin);

            // Register with SPM
            SPMRegistration<float> spm_reg;
            spm_reg.set_reference_image(ref_aladin);
            spm_reg.add_floating_image(floating_sptr);
            spm_reg.add_floating_image(floating_sptr);
            spm_reg.set_working_folder(spm_working_folder);
            spm_reg.set_working_folder_file_overwrite(true);
            spm_reg.set_delete_temp_files(false);
            spm_reg.process();
            const std::shared_ptr<const AffineTransformation<float> > spm_tm_sptr = spm_reg.get_transformation_matrix_forward_sptr(1);
            const AffineTransformation<float> spm_inv_tm = spm_tm_sptr->get_inverse();

            // Check tm roughly equals inverse TM of the resampler
            const std::array<float,3> estimated_euler_angles = spm_inv_tm.get_Euler_angles();
            const std::array<float,3> estimated_translations = { spm_inv_tm[0][3], spm_inv_tm[1][3], spm_inv_tm[2][3] };

            const std::array<float,3> input_euler_angles = tm_sptr->get_Euler_angles();
            const std::array<float,3> input_translations = { (*tm_sptr)[0][3], (*tm_sptr)[1][3], (*tm_sptr)[2][3] };

            std::array<float,3> diff_euler_angles, diff_translations;
            for (unsigned i=0; i<3; ++i) {
                diff_euler_angles[i] = 100.f * (input_euler_angles[i] - estimated_euler_angles[i]) / input_euler_angles[i];
                diff_translations[i] = 100.f * (input_translations[i] - estimated_translations[i]) / input_translations[i];
            }

            std::cout << "Input Euler angles:              " << input_euler_angles[0]     << " " << input_euler_angles[1]     << " " << input_euler_angles[2]     << "\n";
            std::cout << "Estimated Euler angles:          " << estimated_euler_angles[0] << " " << estimated_euler_angles[1] << " " << estimated_euler_angles[2] << "\n";
            std::cout << "Percentage diff in Euler angles: " << diff_euler_angles[0]      << " " << diff_euler_angles[1]      << " " << diff_euler_angles[2]      << "\n";
            std::cout << "Input translations:              " << input_translations[0]     << " " << input_translations[1]     << " " << input_translations[2]     << "\n";
            std::cout << "Estimated translations:          " << estimated_translations[0] << " " << estimated_translations[1] << " " << estimated_translations[2] << "\n";
            std::cout << "Percentage diff in translations: " << diff_translations[0]      << " " << diff_translations[1]      << " " << diff_translations[2]      << "\n";

            // Check differences are less than 1%
            for (unsigned i=0; i<3; ++i) {
                if (std::abs(diff_euler_angles[i]) > 1.f)
                    throw std::runtime_error("SPM registration failed (angles).");
                if (std::abs(diff_translations[i]) > 1.f)
                    throw std::runtime_error("SPM registration failed (translations).");
            }

            if (std::dynamic_pointer_cast<const NiftiImageData3D<float> >(spm_reg.get_output_sptr(1))->operator!=(*ref_aladin))
                throw std::runtime_error("SPM registration failed (image difference).");

        }
        {
            // Try to register via filename
            SPMRegistration<float> spm_reg2;
            spm_reg2.set_reference_image_filename(save_nifti_image);
            spm_reg2.add_floating_image_filename(save_nifti_image);
            spm_reg2.add_floating_image_filename(save_nifti_image);
            spm_reg2.set_working_folder(spm_working_folder2);
            spm_reg2.set_working_folder_file_overwrite(true);
            spm_reg2.set_delete_temp_files(false);
            spm_reg2.process();

            for (unsigned i=0; i<2; ++i) {
                spm_reg2.get_output_sptr(i)->write(output_prefix + "spm_out_" + std::to_string(i));
                spm_reg2.get_displacement_field_forward_sptr(i)->write(output_prefix + "spm_disp_fwd_" + std::to_string(i));
                spm_reg2.get_displacement_field_inverse_sptr(i)->write(output_prefix + "spm_disp_inv_" + std::to_string(i));
                spm_reg2.get_deformation_field_forward_sptr(i)->write(output_prefix + "spm_def_fwd_" + std::to_string(i));
                spm_reg2.get_deformation_field_inverse_sptr(i)->write(output_prefix + "spm_def_inv_" + std::to_string(i));
                spm_reg2.get_transformation_matrix_forward_sptr(i)->write(output_prefix + "spm_tm_fwd_" + std::to_string(i));
                spm_reg2.get_transformation_matrix_inverse_sptr(i)->write(output_prefix + "spm_tm_inv_" + std::to_string(i));
            }
        }

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished SPM test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }
#endif

    // Error handling
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
