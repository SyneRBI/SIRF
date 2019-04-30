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
\brief Class for quaternions.
\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include <array>
#include <vector>

namespace sirf {

// Forward declarations
template<class dataType> class AffineTransformation;

/*!
\ingroup Registration
\brief Class for quaternions.
\author Richard Brown
\author CCP PETMR
*/
template<class dataType>
class Quaternion
{
public: 

    /// Constructor
    Quaternion(const dataType in_w, const dataType in_x, const dataType in_y, const dataType in_z) :
        w(in_w), x(in_x), y(in_y), z(in_z)
    {}

    /// Constructor from transformation matrix
    /// Code from here: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    Quaternion(const AffineTransformation<dataType> tm);

    /// Equality operator
    bool operator==(const Quaternion &other) const;

    /// Equality operator
    bool operator!=(const Quaternion &other) const;

    /// Average quaternions
    static Quaternion get_average(const std::vector<Quaternion> &quaternions);

    /// Dot product with another quaternion
    dataType dot(const Quaternion &other) const;

    /// Changes the sign of the quaternion components. This is not the same as the inverse.
    Quaternion inverse_sign_quaternion() const;

    /// Normalise quaternion
    Quaternion normalise() const;

    /// Print quaternion
    void print() const;

    /// Get the data
    std::array<dataType,4> get_data() const;

    /// Set the data
    void set_data(const std::array<dataType,4> &data);

    /// Get as Euler angles
    std::array<dataType,3> get_Euler_angles() const;

private:

    /// Data
    dataType w, x, y, z;

    /// Returns true if the input quaternion is close to the original. This can
    /// be used to check whether or not one of two quaternions which are supposed to
    /// be very similar but has its component signs reversed (q has the same rotation as -q).
    bool is_quaternion_close(const Quaternion &other) const;
};
}
