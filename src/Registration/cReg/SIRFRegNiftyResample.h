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
\brief Resampling class based on nifty resample

If multiple matrices are set, they will be multiplied in the order that they have been added.
i.e., if A was set first, followed by B, then they will be multiplied AB.

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFREGNIFTYRESAMPLE_H_
#define _SIRFREGNIFTYRESAMPLE_H_

#include <nifti1_io.h>
#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>
#include "SIRFRegMisc.h"

class SIRFRegNiftyResample
{
public:

    enum InterpolationType {
        NOTSET           = -1,
        NEARESTNEIGHBOUR =  0,
        LINEAR           =  1,
        CUBICSPLINE      =  3
    };

    /// Constructor
    SIRFRegNiftyResample() { _interpolation_type = NOTSET; }

    /// Destructor
    virtual ~SIRFRegNiftyResample() {}

    /// Set reference .nii image
    void set_reference_image(const nifti_image *reference_image_ptr)
    {
        _reference_image_sptr     = std::make_shared<nifti_image>(*reference_image_ptr);
        _reference_image_filename = "";
    }

    /// Set reference .nii filename
    void set_reference_image_filename(const std::string reference_image_filename)
    {
        _reference_image_filename = reference_image_filename;
        _reference_image_sptr.reset();
    }

    /// Set floating .nii image
    void set_floating_image(const nifti_image *floating_image_ptr)
    {
        _floating_image_sptr     = std::make_shared<nifti_image>(*floating_image_ptr);
        _floating_image_filename = "";
    }

    /// Set floating .nii filename
    void set_floating_image_filename(const std::string floating_image_filename)
    {
        _floating_image_filename = floating_image_filename;
        _floating_image_sptr.reset();
    }

    /// Add transformation matrix
    void add_transformation_matrix(const mat44 *transformation_matrix)
    {
        _transformation_matrices.push_back(std::make_shared<mat44>(*transformation_matrix));
    }

    /// Add transformation matrix filename (4x4 .csv file)
    void add_transformation_matrix_filename(const std::string transformation_matrix_filename)
    {
        std::shared_ptr<mat44> transformation_matrix;
        SIRFRegMisc::open_transformation_matrix(transformation_matrix,transformation_matrix_filename);
        _transformation_matrices.push_back(transformation_matrix);
    }

    /// Clear transformation matrices
    void clear_transformation_matrices()
    {
        // Empty the vector containing the transformation matrices
        _transformation_matrices.resize(0);
    }

    /// Set interpolation type to nearest neighbour
    void set_interpolation_type_to_nearest_neighbour() { _interpolation_type = NEARESTNEIGHBOUR; }

    /// Set interpolation type to linear
    void set_interpolation_type_to_linear() { _interpolation_type = LINEAR; }

    /// Set interpolation type to cubic spline
    void set_interpolation_type_to_cubic_spline() { _interpolation_type = CUBICSPLINE; }

    /// Update
    void update();

    /// Get output
    nifti_image *get_output() const { return _output_image_sptr.get(); }

    /// Save resampled image to file
    void save_resampled_image(const std::string filename) const;

protected:

    /// Check parameters
    virtual void check_parameters();

    /// Set up the transformation matrix
    void set_up_transformation_matrix(mat44 &matrix);

    /// Set up the deformation field image
    void set_up_deformation_field_image(std::shared_ptr<nifti_image> &deformation_field_image_sptr, mat44 matrix);

    /// Set up the output image
    void set_up_output_image();

    boost::filesystem::path              _reference_image_filename;
    std::shared_ptr<nifti_image>         _reference_image_sptr;

    boost::filesystem::path              _floating_image_filename;
    std::shared_ptr<nifti_image>         _floating_image_sptr;

    std::vector<std::shared_ptr<mat44> > _transformation_matrices;

    InterpolationType                    _interpolation_type;

    std::shared_ptr<nifti_image>         _output_image_sptr;
};

#endif
