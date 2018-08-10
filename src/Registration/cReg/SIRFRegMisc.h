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
\brief Generic tools (e.g., opening files)

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFMISC_H_
#define _SIRFMISC_H_

#include <boost/filesystem.hpp>
#include <nifti1_io.h>
#include <iostream>
#include "SIRFImageData.h"

class SIRFImageDataDeformation;

namespace SIRFRegMisc {

    /// Open nifti image
    void open_nifti_image(std::shared_ptr<nifti_image> &image, const boost::filesystem::path filename);

    /// Save nifti image
    void save_nifti_image(nifti_image *image, const std::string filename);

    /// Save nifti image
    void save_nifti_image(std::shared_ptr<nifti_image> image, const std::string filename);

    //! Split multi-component nifti image
    /*! Assume that the multi-component aspect is in dim[5].
        Split the image into that many single-component images
    */
    std::vector<std::shared_ptr<nifti_image> >
        split_multicomponent_nifti_image(std::shared_ptr<nifti_image> input_sptr);

    /// Save a multicomponent nifti image
    void save_multicomponent_nifti_image(std::shared_ptr<nifti_image> input_sptr, const std::string &filename, const bool &split_xyz);

    /// Copy nifti image
    void copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr);

    /// Flip multicomponent image along a given axis
    void flip_multicomponent_image(SIRFImageDataDeformation &im, int dim);

    /// Get cpp from transformation matrix
#if NIFTYREG_VER_1_3
    void get_cpp_from_transformation_matrix(std::shared_ptr<nifti_image> &cpp_sptr, const std::shared_ptr<mat44> &TM_sptr, const std::shared_ptr<nifti_image> &warped_sptr);
#endif
    /// Get def from cpp
    void get_def_from_cpp(SIRFImageDataDeformation &def, const SIRFImageDataDeformation &cpp, const SIRFImageData &ref);

    /// Convert from deformation to displacement field image
    void convert_from_def_to_disp(SIRFImageDataDeformation &im);

    /// Convert from displacement to deformation field image
    void convert_from_disp_to_def(SIRFImageDataDeformation &im);

    /// Multiply image
    void multiply_image(SIRFImageData &output, const SIRFImageData &input, const float &value);

    /// Do nifti images match?
    bool do_nifti_image_match(const SIRFImageData &im1, const SIRFImageData &im2);

    /// Do nifti image elements match?
    template<typename T>
    bool do_nifti_image_elements_match(const std::string &name, const T &elem1, const T &elem2)
    {
        if(fabs(elem1-elem2) < 1.e-7F)
            return true;
        std::cout << "mismatch in " << name << " , (values: " <<  elem1 << " and " << elem2 << ")\n";
        return false;
    }

    /// Do nifti image elements match?
    bool do_nifti_image_elements_match(const std::string &name, const mat44 &elem1, const mat44 &elem2);

    /// Dump info of nifti image
    void dump_nifti_info(const std::string &im_filename);

    /// Dump info of nifti image
    void dump_nifti_info(const SIRFImageData &im);

    /// Dump info of multiple nifti images
    void dump_nifti_info(const std::vector<SIRFImageData> &ims);

    /// Save transformation matrix to file
    void save_transformation_matrix(const std::shared_ptr<mat44> &transformation_matrix_sptr, const std::string &filename);

    /// Read transformation matrix from file
    void open_transformation_matrix(std::shared_ptr<mat44> &transformation_matrix_sptr, const std::string &filename);

    /// Print mat44
    void print_mat44(const mat44 &mat);

    /// Print multiple mat44
    void print_mat44(const std::vector<mat44> &mats);

    /// Do mat44 match?
    bool do_mat44_match( const mat44& mat1, const mat44& mat2 );

    /// Mat44 multiplier
    mat44 multiply_mat44(const mat44 &x, const mat44 &y);
}

#endif
