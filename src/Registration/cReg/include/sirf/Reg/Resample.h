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
#include "sirf/Reg/Transformation.h"
#include "sirf/iUtilities/iutilities.h"

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

    /// Set reference image. This is the image that would be the reference if you were doing a forward transformation.
    virtual void set_reference_image(const std::shared_ptr<const ImageData> reference_image_sptr);

    /// Set floating image. This is the image that would be the floating if you were doing a forward transformation.
    virtual void set_floating_image(const std::shared_ptr<const ImageData> floating_image_sptr);

    /// Add transformation
    virtual void add_transformation(const std::shared_ptr<const Transformation<dataType> > transformation_sptr);

    /// Set interpolation type (0=nearest neighbour, 1=linear, 3=cubic, 4=sinc)
    virtual void set_interpolation_type(const enum InterpolationType type);

    /// Set interpolation type to nearest neighbour
    void set_interpolation_type_to_nearest_neighbour() { set_interpolation_type(NEARESTNEIGHBOUR); }

    /// Set interpolation type to linear
    void set_interpolation_type_to_linear() { set_interpolation_type(LINEAR); }

    /// Set interpolation type to cubic spline
    void set_interpolation_type_to_cubic_spline() { set_interpolation_type(CUBICSPLINE); }

    /// Set interpolation type to sinc
    void set_interpolation_type_to_sinc() { set_interpolation_type(SINC); }

    /// Get interpolation type
    const InterpolationType get_interpolation_type() const { return _interpolation_type; }

    /// Set padding value
    void set_padding_value(const float padding_value) { _padding_value = padding_value; }

    /// Process - will call forward
    DEPRECATED virtual void process() = 0;

    /// Get output
    const std::shared_ptr<const ImageData> get_output_sptr() const { return _output_image_sptr; }

    /// Do the forward transformation
    virtual std::shared_ptr<ImageData> forward(const std::shared_ptr<const ImageData> input_sptr) = 0;

    /// Do the forward transformation
    virtual void forward(std::shared_ptr<ImageData> output_sptr, const std::shared_ptr<const ImageData> input_sptr) = 0;

    /// Do the adjoint transformation
    virtual std::shared_ptr<ImageData> adjoint(const std::shared_ptr<const ImageData> input_sptr) = 0;

    /// Do the adjoint transformation
    virtual void adjoint(std::shared_ptr<ImageData> output_sptr, const std::shared_ptr<const ImageData> input_sptr) = 0;

    /// Backward. Alias for Adjoint
    virtual std::shared_ptr<ImageData> backward(const std::shared_ptr<const ImageData> input_sptr);

    /// Backward. Alias for Adjoint
    virtual void backward(std::shared_ptr<ImageData> output_sptr, const std::shared_ptr<const ImageData> input_sptr);

protected:

    /// Set up
    virtual void set_up() = 0;

    /// Set up forward
    virtual void set_up_forward() = 0;

    /// Set up adjoint
    virtual void set_up_adjoint() = 0;

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

    /// Padding value
    float _padding_value = 0;
    bool _need_to_set_up = true;
    bool _need_to_set_up_forward = true;
    bool _need_to_set_up_adjoint = true;
};
}
