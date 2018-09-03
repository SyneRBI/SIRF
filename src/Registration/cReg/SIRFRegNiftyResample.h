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

If multiple transformations are set, they will be used in the order that they have been added.
i.e., Trans3(Trans2(Trans1(x))).

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
#include "SIRFImageData.h"
#include "SIRFImageDataDeformation.h"
#include "SIRFRegTransformation.h"

/// Wrapper around NiftyReg's resample class
class SIRFRegNiftyResample
{
public:

    /// Constructor
    SIRFRegNiftyResample() { _interpolation_type = NOTSET; }

    /// Destructor
    virtual ~SIRFRegNiftyResample() {}

    /// Set reference image
    void set_reference_image(const SIRFImageData &reference_image)
    {
        _reference_image = reference_image;
    }

    /// Set floating image
    void set_floating_image(const SIRFImageData &floating_image)
    {
        _floating_image = floating_image;
    }

    /// Add affine transformation
    void add_transformation_affine(const SIRFRegTransformationAffine &affine);

    /// Add displacement transformation
    void add_transformation_disp(const SIRFRegTransformationDisplacement &disp);

    /// Add deformation transformation
    void add_transformation_def(const SIRFRegTransformationDeformation &def);

    /// Set interpolation type (0=nearest neighbour, 1=linear, 3=cubic, 4=sinc)
    void set_interpolation_type(const int type)
    {
        if      (type == 0) _interpolation_type = NEARESTNEIGHBOUR;
        else if (type == 1) _interpolation_type = LINEAR;
        else if (type == 3) _interpolation_type = CUBICSPLINE;
        else if (type == 4) _interpolation_type = SINC;
        else
            throw std::runtime_error("Invalid interpolation type");
    }

    /// Set interpolation type to nearest neighbour
    void set_interpolation_type_to_nearest_neighbour() { _interpolation_type = NEARESTNEIGHBOUR; }

    /// Set interpolation type to linear
    void set_interpolation_type_to_linear() { _interpolation_type = LINEAR; }

    /// Set interpolation type to cubic spline
    void set_interpolation_type_to_cubic_spline() { _interpolation_type = CUBICSPLINE; }

    /// Set interpolation type to sinc
    void set_interpolation_type_to_sinc() { _interpolation_type = SINC; }

    /// Update
    void update();

    /// Get output
    const SIRFImageData &get_output() const { return _output_image; }

    /// Save resampled image to file
    void save_resampled_image(const std::string filename) const;

protected:

    /// Interpolation type
    enum InterpolationType {
        NOTSET           = -1,
        NEARESTNEIGHBOUR =  0,
        LINEAR           =  1,
        CUBICSPLINE      =  3,
        SINC             =  4
    };

    /// Check parameters
    virtual void check_parameters();

    /// Set up the transformation matrix
    void set_up_transformation_matrix(mat44 &matrix);

    /// Set up the output image
    void set_up_output_image();

    /// Reference image
    SIRFImageData            _reference_image;
    /// Floating image
    SIRFImageData            _floating_image;

    /// Transformations (could be mixture of affine, displacements, deformations).
    std::vector<std::shared_ptr<SIRFRegTransformation> > _transformations;

    /// Interpolation type
    InterpolationType        _interpolation_type;

    /// Output image
    SIRFImageData            _output_image;
};

#endif
