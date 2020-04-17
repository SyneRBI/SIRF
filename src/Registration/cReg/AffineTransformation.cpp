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

#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/Quaternion.h"
#include <_reg_globalTrans.h>
#include <iomanip>
#include <boost/filesystem.hpp>

using namespace sirf;

static void check_folder_exists_if_not_create(const std::string &path)
{
    // If the folder doesn't exist, create it
    boost::filesystem::path path_boost(path);
    if (!boost::filesystem::exists(path_boost.parent_path())) {
        if (path_boost.parent_path().string() != "") {
            std::cout << "\n\tCreating folder: \"" << path_boost.parent_path().string() << "\"\n" << std::flush;
            boost::filesystem::create_directory(path_boost.parent_path());
        }
    }
}

template<class dataType>
AffineTransformation<dataType>::AffineTransformation(const dataType tm[4][4])
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = tm[i][j];
}

template<class dataType>
AffineTransformation<dataType>::AffineTransformation(const AffineTransformation& to_copy)
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = to_copy._tm[i][j];
}

template<class dataType>
AffineTransformation<dataType>& AffineTransformation<dataType>::operator=(const AffineTransformation& to_copy)
{
    // check for self-assignment
    if (this != &to_copy) {
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                _tm[i][j] = to_copy._tm[i][j];
    }
    return *this;
}

template<class dataType>
void AffineTransformation<dataType>::print(const std::vector<AffineTransformation> &mats)
{
    for(int i=0;i<4;i++) {
        std::cout << "\t" << std::left << std::setw(19) << "";
        for(unsigned j=0;j<mats.size();j++) {
            std::ostringstream ss;
            ss << "[" <<
                  std::setprecision(3) << mats[j][i][0] << "," <<
                  std::setprecision(3) << mats[j][i][1] << "," <<
                  std::setprecision(3) << mats[j][i][2] << "," <<
                  std::setprecision(3) << mats[j][i][3] << "] ";
            std::cout << std::left << std::setw(19) << ss.str();
        }
        std::cout << "\n";
    }
}

template<class dataType>
AffineTransformation<dataType>::AffineTransformation()
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = i == j ? 1.F : 0.F;
}

template<class dataType>
AffineTransformation<dataType>::AffineTransformation(const std::string &filename)
{
#ifndef NDEBUG
    std::cout << "\n\nReading transformation matrix from file...\n\n";
#endif
    try{
        std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
        // If the file is greater than 10KB it's exxxtremely unlikely it's a 4x4 TM
        if (in.tellg() > 10000)
            throw std::runtime_error("File over 10KB, assuming it's not a TM");

        std::ifstream file(filename);

        // Check opening worked
        if (!file.good())
            throw std::runtime_error("Cannot find the file: " + filename + ".");

        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                file >> _tm[i][j];
    }
    catch (...) {
        throw std::runtime_error("Failed reading affine matrix (" + filename + ").\n\t");
    }
#ifndef NDEBUG
    std::cout << "\n\nSuccessfully read transformation matrix from file:\n";
    this->print();
#endif
}

template<class dataType>
AffineTransformation<dataType>::AffineTransformation(const mat44 &tm)
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            this->_tm[i][j] = dataType(tm.m[i][j]);
}

template<class dataType>
AffineTransformation<dataType>::AffineTransformation(const std::array<dataType,3> &trans, const Quaternion<dataType> &quat)
{
    // Set bottom row to [0,0,0,1]
    for (unsigned i=0; i<3; ++i)
        _tm[3][i] = 0.f;
    _tm[3][3] = 1.f;

    // Translations are easy
    for (unsigned i=0; i<3; ++i)
        _tm[i][3] = trans[i];

    // Convert quaternions back to rotation matrix
    std::array<dataType,4> quat_data = quat.get_data();
    const dataType a(quat_data[0]);
    const dataType b(quat_data[1]);
    const dataType c(quat_data[2]);
    const dataType d(quat_data[3]);

    _tm[0][0] = dataType(2*pow(a,2) -1 + 2*pow(b,2));
    _tm[0][1] = 2*b*c-2*a*d;
    _tm[0][2] = 2*b*d+2*a*c;

    _tm[1][0] = 2*b*c+2*a*d;
    _tm[1][1] = dataType(2*pow(a,2) -1 + 2*pow(c,2));
    _tm[1][2] = 2*c*d-2*a*b;

    _tm[2][0] = 2*b*d-2*a*c;
    _tm[2][1] = 2*c*d+2*a*b;
    _tm[2][2] = dataType(2*pow(a,2) -1 + 2*pow(d,2));
}

template<class dataType>
AffineTransformation<dataType>::AffineTransformation(const std::array<dataType,3> &trans, const std::array<dataType,3> &euler, const bool degrees)
{
    // Set bottom row to [0,0,0,1]
    for (unsigned i=0; i<3; ++i)
        _tm[3][i] = 0.f;
    _tm[3][3] = 1.f;

    // Set translations
    for (unsigned i=0; i<3; ++i)
        _tm[i][3] = trans[i];

    std::array<dataType,3> euler_rad = euler;

    // if in degrees, convert to radian
    if (degrees)
        for (unsigned i=0; i<3; ++i)
            euler_rad[i] *= dataType(M_PI) / 180.f;

    // Convert euler angles to rotation matrix
    float cu = cos(euler_rad[0]);
    float su = sin(euler_rad[0]);
    float cv = cos(euler_rad[1]);
    float sv = sin(euler_rad[1]);
    float cw = cos(euler_rad[2]);
    float sw = sin(euler_rad[2]);

    _tm[0][0] = cv*cw;
    _tm[0][1] = su*sv*cw - cu*sw;
    _tm[0][2] = su*sw + cu*sv*cw;

    _tm[1][0] = cv*sw;
    _tm[1][1] = cu*cw + su*sv*sw;
    _tm[1][2] = cu*sv*sw - su*cw;

    _tm[2][0] = -sv;
    _tm[2][1] = su*cv;
    _tm[2][2] = cu*cv;
}

template<class dataType>
bool AffineTransformation<dataType>::operator==(const AffineTransformation &other) const
{
    if (this == &other)
        return true;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            if( fabs(_tm[j][i] - other._tm[j][i]) > 1.e-4F )
                return false;
    return true;
}

/// Equality operator

template<class dataType>
bool AffineTransformation<dataType>::operator!=(const AffineTransformation &other) const
{
    return !(*this == other);
}

/// Multiply matrices
template<class dataType>
AffineTransformation<dataType> AffineTransformation<dataType>::operator* (const AffineTransformation &other) const
{
    // Create result, set to zero (initially identity)
    AffineTransformation res;
    for (int i=0;i<4;i++)
        res[i][i] = 0.F;

    for (int i=0;i<4;i++)
        for (int j=0;j<4;j++)
            for (int k=0;k<4;k++)
                res[i][j] += (*this)[i][k] * other[k][j];

    return res;
}

template<class dataType>
mat44 AffineTransformation<dataType>::get_as_mat44() const
{
    mat44 tm;
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            tm.m[i][j] = float(_tm[i][j]);
    return tm;
}

template<class dataType>
NiftiImageData3DDeformation<dataType> AffineTransformation<dataType>::get_as_deformation_field(const NiftiImageData<dataType> &ref) const
{
    NiftiImageData3DDeformation<dataType> def;
    def.create_from_3D_image(ref);
    mat44 mat = this->get_as_mat44();
    reg_affine_getDeformationField(&mat, def.get_raw_nifti_sptr().get());
    def.get_raw_nifti_sptr()->intent_p1 = DEF_FIELD;
    return def;
}

template<class dataType>
AffineTransformation<dataType> AffineTransformation<dataType>::deep_copy() const
{
    AffineTransformation temp;
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            temp[i][j] = _tm[i][j];
    return temp;
}

/// Save transformation matrix to file
template<class dataType>
void AffineTransformation<dataType>::write(const std::string &filename) const
{
    // Check that input isn't blank
    if (filename.empty())
        throw std::runtime_error("Error, cannot write transformation matrix to file because filename is blank");

    std::string filename_w_ext = filename;
    boost::filesystem::path filename_boost(filename_w_ext);
    if (!filename_boost.has_extension())
        filename_w_ext += ".txt";

    // If the folder doesn't exist, create it
    check_folder_exists_if_not_create(filename_w_ext);

    FILE *file;
    file=fopen(filename_w_ext.c_str(), "w");
    for(int i=0; i<4; ++i)
        fprintf(file, "%e %e %e %e\n", _tm[i][0], _tm[i][1], _tm[i][2], _tm[i][3]);
    fclose(file);
}

template<class dataType>
dataType AffineTransformation<dataType>::get_determinant() const
{
    return  _tm[0][0]*(_tm[1][1]*_tm[2][2] - _tm[2][1]*_tm[1][2]) -
            _tm[0][1]*(_tm[1][0]*_tm[2][2] - _tm[2][0]*_tm[1][2]) +
            _tm[0][2]*(_tm[1][0]*_tm[2][1] - _tm[2][0]*_tm[1][1]);
}

template<class dataType>
void AffineTransformation<dataType>::print() const
{
    AffineTransformation<dataType>::print({*this});
}

template<class dataType>
AffineTransformation<dataType> AffineTransformation<dataType>::get_inverse() const
{
    mat44 res = nifti_mat44_inverse(this->get_as_mat44());
    return AffineTransformation(res.m);
}

template<class dataType>
const std::array<dataType,3> AffineTransformation<dataType>::get_Euler_angles() const
{
    if (!this->is_rigid())
        throw std::runtime_error("Transformation matrix needs to be rigid in order for Euler angles to be calculated.");

    float sy = sqrt(_tm[0][0] * _tm[0][0] +  _tm[1][0] * _tm[1][0] );
    bool singular = sy < 1e-6F;

    dataType x, y, z;
    if (!singular) {
        x = atan2(_tm[2][1], _tm[2][2]);
        y = atan2(-_tm[2][0], sy);
        z = atan2(_tm[1][0], _tm[0][0]);
    }
    else {
        x = atan2(-_tm[1][2], _tm[1][1]);
        y = atan2(-_tm[2][0], sy);
        z = 0;
    }
    return std::array<dataType,3>{x, y, z};
}

template<class dataType>
Quaternion<dataType> AffineTransformation<dataType>::get_quaternion() const
{
    return Quaternion<dataType>(*this);
}

template<class dataType>
bool AffineTransformation<dataType>::is_rigid() const
{
    return std::abs(dataType(std::pow(this->get_determinant(),2)) - 1.F) < 1.e-4F;
}

template<class dataType>
AffineTransformation<dataType> AffineTransformation<dataType>::get_average(const std::vector<AffineTransformation<dataType> > &mats)
{
    // Array for translations and vector of quaternions
    std::array<dataType,3> avg_trans{0.F,0.F,0.F};
    std::vector<Quaternion<dataType> > quaternions;

    // loop over all matrices
    for (size_t i=0; i<mats.size(); ++i) {
        // sum translations
        for (unsigned j=0; j<3; ++j)
            avg_trans[j] += mats[i][j][3];

        // For quaternions, need to extract them from TM
        quaternions.push_back(Quaternion<dataType>(mats[i]));
    }
    // For translations, average by dividing by number of TMs
    for (unsigned i=0; i<3; ++i)
        avg_trans[i] /= dataType(mats.size());

    // For quaternions, slightly more complicated
    Quaternion<dataType> avg_quaternion = Quaternion<dataType>::get_average(quaternions);

    // Return the average
    return AffineTransformation<dataType>(avg_trans,avg_quaternion);
}

namespace sirf {
template class AffineTransformation<float>;
}
