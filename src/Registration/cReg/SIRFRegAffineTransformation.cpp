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

#include "SIRFRegAffineTransformation.h"
#include "NiftiImageData3DDeformation.h"
#include <_reg_globalTrans.h>
#include <iomanip>

using namespace sirf;

SIRFRegAffineTransformation::SIRFRegAffineTransformation(const float tm[4][4])
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = tm[i][j];
}

SIRFRegAffineTransformation::SIRFRegAffineTransformation(const SIRFRegAffineTransformation& to_copy)
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = to_copy._tm[i][j];
}

SIRFRegAffineTransformation& SIRFRegAffineTransformation::operator=(const SIRFRegAffineTransformation& to_copy)
{
    // check for self-assignment
    if (this != &to_copy) {
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                _tm[i][j] = to_copy._tm[i][j];
    }
    return *this;
}

void SIRFRegAffineTransformation::print(const std::vector<SIRFRegAffineTransformation> &mats)
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

SIRFRegAffineTransformation SIRFRegAffineTransformation::get_identity()
{
    SIRFRegAffineTransformation res;
    for (int i=0; i<4; ++i)
        res[i][i] = 1.F;
    return res;
}

SIRFRegAffineTransformation::SIRFRegAffineTransformation()
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            _tm[i][j] = 0.F;
}

SIRFRegAffineTransformation::SIRFRegAffineTransformation(const std::string &filename)
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

bool SIRFRegAffineTransformation::operator==(const SIRFRegAffineTransformation &other) const
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
bool SIRFRegAffineTransformation::operator!=(const SIRFRegAffineTransformation &other) const
{
    return !(*this == other);
}

/// Multiply matrices
SIRFRegAffineTransformation SIRFRegAffineTransformation::operator* (const SIRFRegAffineTransformation &other) const
{
    // Print info
    std::cout << "\nMultiplying two matrices...\n";
    std::cout << "Matrix 1:\n";
    this->print();
    std::cout << "Matrix 2:\n";
    other.print();

    // Create result, set to zero
    SIRFRegAffineTransformation res;
    for (int i=0;i<4;i++)
        for (int j=0;j<4;j++)
            for (int k=0;k<4;k++)
                res[i][j] += (*this)[i][k] * other[k][j];

    std::cout << "Result:\n";
    res.print();

    return res;
}

mat44 SIRFRegAffineTransformation::get_as_mat44() const
{
    mat44 tm;
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            tm.m[i][j] = _tm[i][j];
    return tm;
}

NiftiImageData3DDeformation SIRFRegAffineTransformation::get_as_deformation_field(const NiftiImageData3D &ref) const
{
    NiftiImageData3DDeformation def;
    def.create_from_3D_image(ref);
    mat44 mat = this->get_as_mat44();
    reg_affine_getDeformationField(&mat, def.get_raw_nifti_sptr().get());
    def.get_raw_nifti_sptr()->intent_p1 = DEF_FIELD;
    return def;
}

SIRFRegAffineTransformation SIRFRegAffineTransformation::deep_copy() const
{
    SIRFRegAffineTransformation temp;
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            temp[i][j] = _tm[i][j];
    return temp;
}

/// Save transformation matrix to file
void SIRFRegAffineTransformation::save_to_file(const std::string &filename) const
{
    // Check that input isn't blank
    if (filename.empty())
        throw std::runtime_error("Error, cannot write transformation matrix to file because filename is blank");

    // If the folder doesn't exist, create it
    SIRFRegMisc::check_folder_exists(filename);

    FILE *file;
    file=fopen(filename.c_str(), "w");
    for(int i=0; i<4; ++i)
        fprintf(file, "%.7F %.7F %.7F %.7F\n", _tm[i][0], _tm[i][1], _tm[i][2], _tm[i][3]);
    fclose(file);
}

float SIRFRegAffineTransformation::get_determinant() const
{
    return  _tm[0][0]*(_tm[1][1]*_tm[2][2] - _tm[2][1]*_tm[1][2]) -
            _tm[0][1]*(_tm[1][0]*_tm[2][2] - _tm[2][0]*_tm[1][2]) +
            _tm[0][2]*(_tm[1][0]*_tm[2][1] - _tm[2][0]*_tm[1][1]);
}

void SIRFRegAffineTransformation::print() const
{
    SIRFRegAffineTransformation::print({*this});
}

SIRFRegAffineTransformation SIRFRegAffineTransformation::get_inverse() const
{
    mat44 res = nifti_mat44_inverse(this->get_as_mat44());
    return SIRFRegAffineTransformation(res.m);
}
