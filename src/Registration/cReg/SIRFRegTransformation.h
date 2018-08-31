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

#ifndef _SIRFREGTRANSFORMATION_H
#define _SIRFREGTRANSFORMATION_H

#include "SIRFImageDataDeformation.h"
#include "SIRFRegMisc.h"

/// Abstract base class for SIRFReg transformations
class SIRFRegTransformation
{
public:

    /// Constructor
    SIRFRegTransformation() {}

    /// Destructor
    virtual ~SIRFRegTransformation() {}

    /// Get as deformation field
    virtual SIRFImageDataDeformation get_as_deformation_field(const SIRFImageData &) = 0;
};

/// Class for SIRFReg transformations with an affine transformation
class SIRFRegTransformationAffine : public SIRFRegTransformation
{
public:

    /// Default constructor
    SIRFRegTransformationAffine() {}

    /// Constructor
    SIRFRegTransformationAffine(const mat44 &tm) { _tm = tm; }

    /// Construct from file
    SIRFRegTransformationAffine(const std::string &filename)
    {
        SIRFRegMisc::open_transformation_matrix(_tm,filename);
    }

    /// Destructor
    virtual ~SIRFRegTransformationAffine() {}

    /// Get as deformation field
    virtual SIRFImageDataDeformation get_as_deformation_field(const SIRFImageData &ref);

    /// Deep copy
    virtual SIRFRegTransformationAffine deep_copy() const;

protected:
    mat44 _tm;
};

/// Class for SIRFReg transformations with a displacement field image
class SIRFRegTransformationDisplacement : public SIRFRegTransformation
{
public:

    /// Default constructor
    SIRFRegTransformationDisplacement() {}

    /// Constructor
    SIRFRegTransformationDisplacement(const SIRFImageDataDeformation &disp) { _disp = disp.deep_copy(); }

    /// Construct from file
    SIRFRegTransformationDisplacement(const std::string &filename) { _disp = SIRFImageDataDeformation(filename).deep_copy(); }

    /// Destructor
    virtual ~SIRFRegTransformationDisplacement() {}

    /// Get as deformation field
    virtual SIRFImageDataDeformation get_as_deformation_field(const SIRFImageData &);

    /// Deep copy
    virtual SIRFRegTransformationDisplacement deep_copy() const;


protected:
    SIRFImageDataDeformation _disp;
};

/// Class for SIRFReg transformations with a deformation field image
class SIRFRegTransformationDeformation : public SIRFRegTransformation
{
public:

    /// Default constructor
    SIRFRegTransformationDeformation() {}

    /// Constructor
    SIRFRegTransformationDeformation(const SIRFImageDataDeformation &def) { _def = def.deep_copy(); }

    /// Construct from file
    SIRFRegTransformationDeformation(const std::string &filename) { _def = SIRFImageDataDeformation(filename).deep_copy(); }

    /// Destructor
    virtual ~SIRFRegTransformationDeformation() {}

    /// Get as deformation field
    virtual SIRFImageDataDeformation get_as_deformation_field(const SIRFImageData &);

    /// Deep copy
    virtual SIRFRegTransformationDeformation deep_copy() const;


protected:
    SIRFImageDataDeformation _def;
};


#endif
