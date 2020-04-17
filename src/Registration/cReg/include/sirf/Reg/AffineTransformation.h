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
\brief Class for affine transformations.
\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include <vector>
#include <nifti1_io.h>
#include "sirf/Reg/Transformation.h"
#include <array>

namespace sirf {

// Forward declarations
template<class dataType> class Quaternion;

/*!
\ingroup Registration
\brief Class for affine transformations.
\author Richard Brown
\author CCP PETMR
*/
template<class dataType>
class AffineTransformation : public Transformation<dataType>
{
public:

    /// Print multiple AffineTransformation
    static void print(const std::vector<AffineTransformation<dataType> > &mats);

    /// Default constructor - identity matrix
    AffineTransformation();

    /// Constructor
    AffineTransformation(const dataType tm[4][4]);

    /// Construct from file
    AffineTransformation(const std::string &filename);

    /// Construct from mat44
    AffineTransformation(const mat44 &tm);

    /// Construct from translation and quaternion
    /// Code from here: https://uk.mathworks.com/help/robotics/ref/quaternion.rotmat.html
    AffineTransformation(const std::array<dataType,3> &trans, const Quaternion<dataType> &quat);

    /// Construct from translation and euler angles (XYZ order)
    AffineTransformation(const std::array<dataType,3> &trans, const std::array<dataType,3> &euler, const bool degrees = true);

    /// Copy constructor
    AffineTransformation(const AffineTransformation& to_copy);

    /// Assignment
    AffineTransformation& operator=(const AffineTransformation& to_copy);

    /// Equality operator
    bool operator==(const AffineTransformation &other) const;

    /// Equality operator
    bool operator!=(const AffineTransformation &other) const;

    /// Multiplication operator
    AffineTransformation operator* (const AffineTransformation &other) const;

    /// Overload [] operator (const)
    dataType const *operator [](unsigned i) const { return _tm[i]; }

    /// Overload [] operator
    dataType *operator [](unsigned i) { return _tm[i]; }

    /// Get raw mat44
    mat44 get_as_mat44() const;

    /// Destructor
    virtual ~AffineTransformation() {}

    /// Get as deformation field
    virtual NiftiImageData3DDeformation<dataType> get_as_deformation_field(const NiftiImageData<dataType> &ref) const;

    /// Deep copy
    virtual AffineTransformation deep_copy() const;

    /// Save to file
    virtual void write(const std::string &filename) const;

    /// Get determinant
    dataType get_determinant() const;

    /// Print
    void print() const;

    /// Get inverse
    AffineTransformation get_inverse() const;

    /// Get Euler angles (XYZ)
    const std::array<dataType,3> get_Euler_angles() const;

    /// Get quaternion
    Quaternion<dataType> get_quaternion() const;

    /// Is rigid? If so, determinant will be +/- 1.
    bool is_rigid() const;

    /// Average transformation matrices (using quaternions for rotation component)
    static AffineTransformation get_average(const std::vector<AffineTransformation<dataType> > &mats);

protected:
    dataType _tm[4][4];
};
}
