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
\brief Classes for SIRFReg transformations.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/cReg/AffineTransformation.h"
#include "sirf/cReg/NiftiImageData3DDeformation.h"
#include <_reg_globalTrans.h>
#include <iomanip>

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
SIRFRegAffineTransformation<dataType>::SIRFRegAffineTransformation(const dataType tm[4][4])
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = tm[i][j];
}

template<class dataType>
SIRFRegAffineTransformation<dataType>::SIRFRegAffineTransformation(const SIRFRegAffineTransformation& to_copy)
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = to_copy._tm[i][j];
}

template<class dataType>
SIRFRegAffineTransformation<dataType>& SIRFRegAffineTransformation<dataType>::operator=(const SIRFRegAffineTransformation& to_copy)
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
void SIRFRegAffineTransformation<dataType>::print(const std::vector<SIRFRegAffineTransformation> &mats)
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
SIRFRegAffineTransformation<dataType>::SIRFRegAffineTransformation()
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = i == j ? 1.F : 0.F;
}

template<class dataType>
SIRFRegAffineTransformation<dataType>::SIRFRegAffineTransformation(const std::string &filename)
{
    // Check that the file exists
    if (!boost::filesystem::exists(filename))
        throw std::runtime_error("Cannot find the file: " + filename + ".");

    std::cout << "\n\nReading transformation matrix from file...\n\n";

    try{
        std::ifstream file(filename);
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                file >> _tm[i][j];
    }
    catch (...) {
        throw std::runtime_error("Failed reading affine matrix (" + filename + ").\n\t");
    }

    std::cout << "\n\nSuccessfully read transformation matrix from file:\n";

    this->print();
}

template<class dataType>
SIRFRegAffineTransformation<dataType>::SIRFRegAffineTransformation(const mat44 &tm)
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            this->_tm[i][j] = dataType(tm.m[i][j]);
}

template<class dataType>
bool SIRFRegAffineTransformation<dataType>::operator==(const SIRFRegAffineTransformation &other) const
{
    if (this == &other)
        return true;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            if( fabs(_tm[j][i] - other._tm[j][i]) > 1.e-5F )
                return false;
    return true;
}

/// Equality operator

template<class dataType>
bool SIRFRegAffineTransformation<dataType>::operator!=(const SIRFRegAffineTransformation &other) const
{
    return !(*this == other);
}

/// Multiply matrices
template<class dataType>
SIRFRegAffineTransformation<dataType> SIRFRegAffineTransformation<dataType>::operator* (const SIRFRegAffineTransformation &other) const
{
    // Create result, set to zero (initially identity)
    SIRFRegAffineTransformation res;    
    for (int i=0;i<4;i++)
        res[i][i] = 0.F;

    for (int i=0;i<4;i++)
        for (int j=0;j<4;j++)
            for (int k=0;k<4;k++)
                res[i][j] += (*this)[i][k] * other[k][j];

#ifndef NDEBUG
    std::cout << "\nMultiplying two matrices...\n";
    std::cout << "Matrix 1:\n";
    this->print();
    std::cout << "Matrix 2:\n";
    other.print();
    std::cout << "Result:\n";
    res.print();
#endif
    return res;
}

template<class dataType>
mat44 SIRFRegAffineTransformation<dataType>::get_as_mat44() const
{
    mat44 tm;
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            tm.m[i][j] = float(_tm[i][j]);
    return tm;
}

template<class dataType>
NiftiImageData3DDeformation<dataType> SIRFRegAffineTransformation<dataType>::get_as_deformation_field(const NiftiImageData3D<dataType> &ref) const
{
    NiftiImageData3DDeformation<dataType> def;
    def.create_from_3D_image(ref);
    mat44 mat = this->get_as_mat44();
    reg_affine_getDeformationField(&mat, def.get_raw_nifti_sptr().get());
    def.get_raw_nifti_sptr()->intent_p1 = DEF_FIELD;
    return def;
}

template<class dataType>
SIRFRegAffineTransformation<dataType> SIRFRegAffineTransformation<dataType>::deep_copy() const
{
    SIRFRegAffineTransformation temp;
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            temp[i][j] = _tm[i][j];
    return temp;
}

/// Save transformation matrix to file
template<class dataType>
void SIRFRegAffineTransformation<dataType>::write(const std::string &filename) const
{
    // Check that input isn't blank
    if (filename.empty())
        throw std::runtime_error("Error, cannot write transformation matrix to file because filename is blank");

    // If the folder doesn't exist, create it
    check_folder_exists_if_not_create(filename);

    FILE *file;
    file=fopen(filename.c_str(), "w");
    for(int i=0; i<4; ++i)
        fprintf(file, "%e %e %e %e\n", _tm[i][0], _tm[i][1], _tm[i][2], _tm[i][3]);
    fclose(file);
}

template<class dataType>
dataType SIRFRegAffineTransformation<dataType>::get_determinant() const
{
    return  _tm[0][0]*(_tm[1][1]*_tm[2][2] - _tm[2][1]*_tm[1][2]) -
            _tm[0][1]*(_tm[1][0]*_tm[2][2] - _tm[2][0]*_tm[1][2]) +
            _tm[0][2]*(_tm[1][0]*_tm[2][1] - _tm[2][0]*_tm[1][1]);
}

template<class dataType>
void SIRFRegAffineTransformation<dataType>::print() const
{
    SIRFRegAffineTransformation<dataType>::print({*this});
}

template<class dataType>
SIRFRegAffineTransformation<dataType> SIRFRegAffineTransformation<dataType>::get_inverse() const
{
    mat44 res = nifti_mat44_inverse(this->get_as_mat44());
    return SIRFRegAffineTransformation(res.m);
}

namespace sirf {
template class SIRFRegAffineTransformation<float>;
}
