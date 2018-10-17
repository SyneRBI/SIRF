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

#ifndef _SIRFREGMAT44_H_
#define _SIRFREGMAT44_H_

#include "SIRFRegTransformation.h"

namespace sirf {
/// Class for SIRFReg transformations with an affine transformation(mat44)
class SIRFRegMat44 : public SIRFRegTransformation
{
public:
    /// Print multiple SIRFRegMat44
    static void print(const std::vector<sirf::SIRFRegMat44> &mats);

    /// Get as identity matrix
    static SIRFRegMat44 get_identity();

    /// Default constructor
    SIRFRegMat44();

    /// Constructor
    SIRFRegMat44(const mat44 &tm) { _tm = tm; }

    /// Construct from file
    SIRFRegMat44(const std::string &filename);

    /// Equality operator
    bool operator==(const SIRFRegMat44 &other) const;

    /// Equality operator
    bool operator!=(const SIRFRegMat44 &other) const;

    /// Multiplication operator
    SIRFRegMat44 operator* (const SIRFRegMat44 &other) const;

    /// Overload [] operator (const)
    float const *operator [](int i) const { return _tm.m[i]; }

    /// Overload [] operator
    float *operator [](int i) { return _tm.m[i]; }

    /// Get raw mat44
    const mat44 &get_raw_mat44() const { return _tm; }

    /// Destructor
    virtual ~SIRFRegMat44() {}

    /// Get as deformation field
    virtual NiftiImageData3DDeformation get_as_deformation_field(const NiftiImageData3D &ref) const;

    /// Deep copy
    virtual SIRFRegMat44 deep_copy() const;

    /// Save to file
    virtual void save_to_file(const std::string &filename) const;

    /// Get determinant
    float get_determinant() const;

    /// Print
    void print() const;

protected:
    mat44 _tm;
};
}

#endif
