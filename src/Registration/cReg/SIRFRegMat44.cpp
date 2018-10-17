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

#include "SIRFRegMat44.h"
#include "NiftiImageData3DDeformation.h"
#include <_reg_ReadWriteMatrix.h>
#include <_reg_globalTrans.h>
#include <iomanip>

using namespace sirf;

void SIRFRegMat44::print(const std::vector<SIRFRegMat44> &mats)
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

SIRFRegMat44 SIRFRegMat44::get_identity()
{
    SIRFRegMat44 res;
    for (int i=0; i<4; ++i)
        res[i][i] = 1.F;
    return res;
}

SIRFRegMat44::SIRFRegMat44()
{
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            (*this)[i][j] = 0.F;
}

SIRFRegMat44::SIRFRegMat44(const std::string &filename)
{
    // Check that the file exists
    if (!boost::filesystem::exists(filename))
        throw std::runtime_error("Cannot find the file: " + filename + ".");

    std::cout << "\n\nReading transformation matrix from file...\n\n";
    reg_tool_ReadAffineFile(&_tm, const_cast<char*>(filename.c_str()));
    std::cout << "\n\nSuccessfully read transformation matrix from file:\n";

    this->print();
}

bool SIRFRegMat44::operator==(const SIRFRegMat44 &other) const
{
    if (this == &other)
        return true;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            if( fabs(_tm.m[j][i] - other._tm.m[j][i]) > 1.e-5F )
                return false;
    return true;
}

/// Equality operator
bool SIRFRegMat44::operator!=(const SIRFRegMat44 &other) const
{
    return !(*this == other);
}

/// Mat44 multiplier
SIRFRegMat44 SIRFRegMat44::operator* (const SIRFRegMat44 &other) const
{
    // Print info
    std::cout << "\nMultiplying two matrices...\n";
    std::cout << "Matrix 1:\n";
    this->print();
    std::cout << "Matrix 2:\n";
    other.print();

    // Create result, set to zero
    SIRFRegMat44 res;
    for (int i=0;i<4;i++)
        for (int j=0;j<4;j++)
            for (int k=0;k<4;k++)
                res[i][j] += (*this)[i][k] * other[k][j];

    std::cout << "Result:\n";
    res.print();

    return res;
}

NiftiImageData3DDeformation SIRFRegMat44::get_as_deformation_field(const NiftiImageData3D &ref) const
{
    NiftiImageData3DDeformation def;
    def.create_from_3D_image(ref);
    mat44 temp = _tm; // Need temp as the following isn't marked const
    reg_affine_getDeformationField(&temp, def.get_raw_nifti_sptr().get());
    def.get_raw_nifti_sptr()->intent_p1 = DEF_FIELD;
    return def;
}

SIRFRegMat44 SIRFRegMat44::deep_copy() const
{
    SIRFRegMat44 temp(_tm);
    return temp;
}

/// Save transformation matrix to file
void SIRFRegMat44::save_to_file(const std::string &filename) const
{
    // Check that input isn't blank
    if (filename.empty())
        throw std::runtime_error("Error, cannot write transformation matrix to file because filename is blank");
    // Need to copy the tm, since the function is not marked const
    mat44 temp = _tm;
    reg_tool_WriteAffineFile(&temp, filename.c_str());
}

float SIRFRegMat44::get_determinant() const
{
    return  _tm.m[0][0]*(_tm.m[1][1]*_tm.m[2][2] - _tm.m[2][1]*_tm.m[1][2]) -
            _tm.m[0][1]*(_tm.m[1][0]*_tm.m[2][2] - _tm.m[2][0]*_tm.m[1][2]) +
            _tm.m[0][2]*(_tm.m[1][0]*_tm.m[2][1] - _tm.m[2][0]*_tm.m[1][1]);
}

void SIRFRegMat44::print() const
{
    SIRFRegMat44::print({*this});
}
