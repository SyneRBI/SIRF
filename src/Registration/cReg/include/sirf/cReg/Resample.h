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
\brief Abstract resampling base class

\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include <vector>
#include <memory>
#include "sirf/cReg/Transformation.h"

namespace sirf {

// Forward declarations
template<class dataType> class Transformation;
class ImageData;

/*!
\file
\ingroup Registration
\brief Abstract resampling base class

If multiple transformations are set, they will be used in the order that they have been added.
i.e., Trans3(Trans2(Trans1(x))).

\author Richard Brown
\author CCP PETMR
*/
template<class dataType>
class Resample
{
public:

    /// Interpolation type
    enum InterpolationType {
        NOTSET           = -1,
        NEARESTNEIGHBOUR =  0,
        LINEAR           =  1,
        CUBICSPLINE      =  3,
        SINC             =  4
    };

    /// Constructor
    Resample() { _interpolation_type = NOTSET; }

    /// Destructor
    virtual ~Resample() {}

    /// Set reference image
    virtual void set_reference_image(const std::shared_ptr<const ImageData> reference_image_sptr) { _reference_image_sptr = reference_image_sptr; }

    /// Set floating image
    virtual void set_floating_image(const std::shared_ptr<const ImageData> floating_image_sptr) { _floating_image_sptr = floating_image_sptr; }

    /// Add transformation
    virtual void add_transformation(const std::shared_ptr<const Transformation<dataType> > transformation_sptr);

    /// Set interpolation type (0=nearest neighbour, 1=linear, 3=cubic, 4=sinc)
    virtual void set_interpolation_type(const enum InterpolationType type)
    {
        _interpolation_type = type;
    }

    /// Set interpolation type to nearest neighbour
    void set_interpolation_type_to_nearest_neighbour() { _interpolation_type = NEARESTNEIGHBOUR; }

    /// Set interpolation type to linear
    void set_interpolation_type_to_linear() { _interpolation_type = LINEAR; }

    /// Set interpolation type to cubic spline
    void set_interpolation_type_to_cubic_spline() { _interpolation_type = CUBICSPLINE; }

    /// Set interpolation type to sinc
    void set_interpolation_type_to_sinc() { _interpolation_type = SINC; }

    /// Process
    virtual void process() = 0;

    /// Get output
    const std::shared_ptr<const ImageData> get_output_sptr() const { return _output_image_sptr; }

protected:

    /// Check parameters
    virtual void check_parameters();

    /// Reference image
    std::shared_ptr<const ImageData> _reference_image_sptr;
    /// Floating image
    std::shared_ptr<const ImageData> _floating_image_sptr;

    /// Transformations (could be mixture of affine, displacements, deformations).
    std::vector<std::shared_ptr<const Transformation<dataType> > > _transformations;

    /// Interpolation type
    InterpolationType  _interpolation_type;

    /// Output image
    std::shared_ptr<ImageData> _output_image_sptr;
};
}
