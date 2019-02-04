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
\brief Test for creation of a weighted-mean image

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include "sirf/cReg/NiftyAladinSym.h"
#include "sirf/cReg/NiftyF3dSym.h"
#include "sirf/cReg/NiftyResample.h"
#include "sirf/cReg/NiftiImageData3D.h"
#include "sirf/cReg/ImageWeightedMean.h"
#include "sirf/cReg/NiftiImageData3DDisplacement.h"
#include "sirf/cReg/AffineTransformation.h"
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
    const std::string flo_aladin_as_unsigned_int                 = output_prefix   + "flo_aladin_as_unsigned_int.nii";
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
        NiftiImageData3D<float>::print_headers({&b, &t});
        for (unsigned i=0; i<b.get_raw_nifti_sptr()->nvox; ++i)
            if (std::abs(b(i) - t(i)) > 1e-4F)
                throw std::runtime_error("NiftiImageData3D constructor from array.");
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
        NiftiImageData3DTensor<float> h(im1, im2, im3);

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

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished NiftiImageData3DDeformation test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    NiftyAladinSym<float> NA;
    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Nifty aladin test...\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        NA.set_reference_image               (            ref_aladin         );
        NA.set_floating_image                (            flo_aladin         );
        NA.set_parameter_file                (      parameter_file_aladin    );
        NA.set_parameter("SetInterpolationToCubic");
        NA.set_parameter("SetLevelsToPerform","1");
        NA.set_parameter("SetMaxIterations","5");
        NA.process();

        // Get outputs
        std::shared_ptr<const NiftiImageData3D<float> >             warped_sptr       = NA.get_output_sptr();
        std::shared_ptr<const AffineTransformation<float> >  TM_forward_sptr   = std::dynamic_pointer_cast<const AffineTransformation<float> > (NA.get_transformation_matrix_forward_sptr());
        std::shared_ptr<const AffineTransformation<float> >  TM_inverse_sptr   = std::dynamic_pointer_cast<const AffineTransformation<float> > (NA.get_transformation_matrix_forward_sptr());
        std::shared_ptr<const NiftiImageData3DDeformation<float> >  def_forward_sptr  = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> > (NA.get_deformation_field_forward_sptr());
        std::shared_ptr<const NiftiImageData3DDeformation<float> >  def_inverse_sptr  = std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> > (NA.get_deformation_field_inverse_sptr());
        std::shared_ptr<const NiftiImageData3DDisplacement<float> > disp_forward_sptr = std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(NA.get_displacement_field_forward_sptr());
        std::shared_ptr<const NiftiImageData3DDisplacement<float> > disp_inverse_sptr = std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(NA.get_displacement_field_inverse_sptr());

        warped_sptr->write    (      aladin_warped    );
        TM_forward_sptr->write(       TM_forward      );
        TM_inverse_sptr->write(       TM_inverse      );
        disp_forward_sptr->write(aladin_disp_forward);
        disp_inverse_sptr->write_split_xyz_components(aladin_disp_inverse);
        def_forward_sptr->write(aladin_def_forward);
        def_inverse_sptr->write_split_xyz_components(aladin_def_inverse);

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

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished Nifty aladin test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    {
        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Starting Nifty f3d test..\n";
        std::cout << "//------------------------------------------------------------------------ //\n";

        NiftyF3dSym<float> NF;
        NF.set_reference_image               (           ref_f3d          );
        NF.set_floating_image                (           flo_f3d          );
        NF.set_parameter_file                (     parameter_file_f3d     );
        NF.set_reference_time_point          (             1              );
        NF.set_floating_time_point           (             1              );
        NF.process();

        // Get outputs
        std::shared_ptr<const NiftiImageData3D<float> >             warped_sptr       = NF.get_output_sptr();
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

        std::cout << "Testing rigid resample...\n";
        NiftyResample<float> nr1;
        nr1.set_reference_image(ref_aladin);
        nr1.set_floating_image(flo_aladin);
        nr1.set_interpolation_type_to_cubic_spline(); // try different interpolations
        nr1.set_interpolation_type(NiftyResample<float>::CUBICSPLINE); // try different interpolations (cubic)
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
        nr2.process();
        nr2.get_output_sptr()->write(nonrigid_resample_disp);

        std::cout << "Testing non-rigid deformation...\n";
        NiftyResample<float> nr3;
        nr3.set_reference_image(ref_aladin);
        nr3.set_floating_image(flo_aladin);
        nr3.set_interpolation_type_to_nearest_neighbour(); // try different interpolations
        nr3.add_transformation(deff);
        nr3.set_interpolation_type_to_linear();
        nr3.process();
        nr3.get_output_sptr()->write(nonrigid_resample_def);

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

        std::cout << "// ----------------------------------------------------------------------- //\n";
        std::cout << "//                  Finished AffineTransformation test.\n";
        std::cout << "//------------------------------------------------------------------------ //\n";
    }

    // Error handling
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
