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
\brief Class for affine transformations.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/Quaternion.h"
#include "sirf/Reg/AffineTransformation.h"
#include <cmath>

using namespace sirf;

template<class dataType>
Quaternion<dataType>::Quaternion(const AffineTransformation<dataType> tm)
{
    dataType trace = tm[0][0] + tm[1][1] + tm[2][2];
    if( trace > 0 ) {
        dataType s = 0.5f / sqrtf(trace+ 1.0f);
        w = 0.25f / s;
        x = ( tm[2][1] - tm[1][2] ) * s;
        y = ( tm[0][2] - tm[2][0] ) * s;
        z = ( tm[1][0] - tm[0][1] ) * s;
    }
    else if ( tm[0][0] > tm[1][1] && tm[0][0] > tm[2][2] ) {
        dataType s = 2.0f * sqrtf( 1.0f + tm[0][0] - tm[1][1] - tm[2][2]);
        w = (tm[2][1] - tm[1][2] ) / s;
        x = 0.25f * s;
        y = (tm[0][1] + tm[1][0] ) / s;
        z = (tm[0][2] + tm[2][0] ) / s;
    }
    else if (tm[1][1] > tm[2][2]) {
        dataType s = 2.0f * sqrtf( 1.0f + tm[1][1] - tm[0][0] - tm[2][2]);
        w = (tm[0][2] - tm[2][0] ) / s;
        x = (tm[0][1] + tm[1][0] ) / s;
        y = 0.25f * s;
        z = (tm[1][2] + tm[2][1] ) / s;
    }
    else {
        dataType s = 2.0f * sqrtf( 1.0f + tm[2][2] - tm[0][0] - tm[1][1] );
        w = (tm[1][0] - tm[0][1] ) / s;
        x = (tm[0][2] + tm[2][0] ) / s;
        y = (tm[1][2] + tm[2][1] ) / s;
        z = 0.25f * s;
    }
    if (w<0)
        *this = this->inverse_sign_quaternion();
}

template<class dataType>
bool Quaternion<dataType>::operator==(const Quaternion &other) const
{
    if (this == &other)
        return true;
    if (std::fabs(w - other.w) < 1e-4f &&
            std::fabs(x - other.x) < 1e-4f &&
            std::fabs(y - other.y) < 1e-4f &&
            std::fabs(z - other.z) < 1e-4f)
        return true;
    return false;
}

template<class dataType>
bool Quaternion<dataType>::operator!=(const Quaternion &other) const
{
    return !(*this == other);
}

template<class dataType>
Quaternion<dataType> Quaternion<dataType>::get_average(const std::vector<Quaternion> &quaternions)
{
    // Check
    if (quaternions.size() == 0)
        throw std::runtime_error("Quaternion<dataType>::average_quaternions: No quaternions given");
    // If only 1 quaternion, average is easy!
    else if (quaternions.size() == 1)
        return quaternions[0];

    // Sum up all quaternions. Start with first, loop over rest
    Quaternion result(quaternions[0]);

    for (unsigned i=1; i<quaternions.size(); ++i) {

        // But first, check whether subsequent quaternions need to be inverted.
        // Because q and -q are the same rotation, but cannot be averaged, we have to make sure they are all the same.
        if (quaternions[i].is_quaternion_close(quaternions[0])) {
            result.w += quaternions[i].w;
            result.x += quaternions[i].x;
            result.y += quaternions[i].y;
            result.z += quaternions[i].z;
        }
        else {
            Quaternion flipped = quaternions[i].inverse_sign_quaternion();
            result.w += flipped.w;
            result.x += flipped.x;
            result.y += flipped.y;
            result.z += flipped.z;
        }
    }

    // Divide by number of quaternions
    result.w /= dataType(quaternions.size());
    result.x /= dataType(quaternions.size());
    result.y /= dataType(quaternions.size());
    result.z /= dataType(quaternions.size());

    return result;
}

template<class dataType>
Quaternion<dataType> Quaternion<dataType>::normalise() const
{
    dataType new_w, new_x, new_y, new_z;
    dataType lengthD = 1.0f / (w*w + x*x + y*y + z*z);
    new_w = w * lengthD;
    new_x = x * lengthD;
    new_y = y * lengthD;
    new_z = z * lengthD;
    return Quaternion<dataType>(new_w, new_x, new_y, new_z);
}

template<class dataType>
Quaternion<dataType> Quaternion<dataType>::inverse_sign_quaternion() const
{
    return Quaternion<dataType>(-w, -x, -y, -z);
}

template<class dataType>
dataType Quaternion<dataType>::dot(const Quaternion<dataType> &other) const
{
    return x*other.x + y*other.y + z*other.z + w*other.w;
}

template<class dataType>
bool Quaternion<dataType>::is_quaternion_close(const Quaternion &other) const
{
    return this->dot(other) >= 0.f;
}

namespace sirf {
template class Quaternion<float>;
}
