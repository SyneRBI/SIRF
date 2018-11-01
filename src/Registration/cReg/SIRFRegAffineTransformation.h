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

#ifndef _SIRFREGAFFINETRANSFORMATION_H_
#define _SIRFREGAFFINETRANSFORMATION_H_

#include "SIRFRegTransformation.h"

namespace sirf {
/// Class for SIRFReg transformations with an affine transformation
class SIRFRegAffineTransformation : public SIRFRegTransformation
{
public:
    /// Print multiple SIRFRegAffineTransformation
    static void print(const std::vector<sirf::SIRFRegAffineTransformation> &mats);

    /// Get as identity matrix
    static SIRFRegAffineTransformation get_identity();

    /// Default constructor
    SIRFRegAffineTransformation();

    /// Constructor
    SIRFRegAffineTransformation(const float tm[4][4]);

    /// Construct from file
    SIRFRegAffineTransformation(const std::string &filename);

    /// Copy constructor
    SIRFRegAffineTransformation(const SIRFRegAffineTransformation& to_copy);

    /// Assignment
    SIRFRegAffineTransformation& operator=(const SIRFRegAffineTransformation& to_copy);

    /// Equality operator
    bool operator==(const SIRFRegAffineTransformation &other) const;

    /// Equality operator
    bool operator!=(const SIRFRegAffineTransformation &other) const;

    /// Multiplication operator
    SIRFRegAffineTransformation operator* (const SIRFRegAffineTransformation &other) const;

    /// Overload [] operator (const)
    float const *operator [](int i) const { return _tm[i]; }

    /// Overload [] operator
    float *operator [](int i) { return _tm[i]; }

    /// Get raw mat44
    mat44 get_as_mat44() const;

    /// Destructor
    virtual ~SIRFRegAffineTransformation() {}

    /// Get as deformation field
    virtual NiftiImageData3DDeformation get_as_deformation_field(const NiftiImageData3D &ref) const;

    /// Deep copy
    virtual SIRFRegAffineTransformation deep_copy() const;

    /// Save to file
    virtual void save_to_file(const std::string &filename) const;

    /// Get determinant
    float get_determinant() const;

    /// Print
    void print() const;

    /// Get inverse
    SIRFRegAffineTransformation get_inverse() const;

protected:
    float _tm[4][4];
};
}

#endif
