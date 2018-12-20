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

#include <vector>
#include <nifti1_io.h>
#include "sirf/cReg/SIRFRegTransformation.h"

namespace sirf {
/// Class for SIRFReg transformations with an affine transformation
template<class dataType>
class SIRFRegAffineTransformation : public SIRFRegTransformation<dataType>
{
public:
    /// Print multiple SIRFRegAffineTransformation
    static void print(const std::vector<sirf::SIRFRegAffineTransformation<dataType> > &mats);

    /// Get as identity matrix
    static SIRFRegAffineTransformation get_identity();

    /// Default constructor
    SIRFRegAffineTransformation();

    /// Constructor
    SIRFRegAffineTransformation(const dataType tm[4][4]);

    /// Construct from file
    SIRFRegAffineTransformation(const std::string &filename);

    /// Construct from mat44
    SIRFRegAffineTransformation(const mat44 &tm);

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
    dataType const *operator [](int i) const { return _tm[i]; }

    /// Overload [] operator
    dataType *operator [](int i) { return _tm[i]; }

    /// Get raw mat44
    mat44 get_as_mat44() const;

    /// Destructor
    virtual ~SIRFRegAffineTransformation() {}

    /// Get as deformation field
    virtual NiftiImageData3DDeformation<dataType> get_as_deformation_field(const NiftiImageData3D<dataType> &ref) const;

    /// Deep copy
    virtual SIRFRegAffineTransformation deep_copy() const;

    /// Save to file
    virtual void write(const std::string &filename) const;

    /// Get determinant
    dataType get_determinant() const;

    /// Print
    void print() const;

    /// Get inverse
    SIRFRegAffineTransformation get_inverse() const;

protected:
    dataType _tm[4][4];
};
}

#endif
