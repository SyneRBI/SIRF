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
\brief Base class for SIRF nifti image data.

\author Richard Brown
\author CCP PETMR
*/

#ifndef _NIFTIIMAGE_H_
#define _NIFTIIMAGE_H_

#include <nifti1_io.h>
#include <string>
#include <memory>
#include <iostream>

namespace sirf {

/// SIRF image data
class NiftiImage
{
public:

    /// Constructor
    NiftiImage() {}

    /// Destructor
    virtual ~NiftiImage() {}

    /// Assignment
    NiftiImage operator=(const NiftiImage& to_copy);

    /// Filename constructor
    NiftiImage(const std::string &filename);

    /// Nifti constructor
    NiftiImage(const nifti_image &image_nifti);

    /// Nifti shared_ptr constructor
    NiftiImage(const std::shared_ptr<nifti_image> image_nifti);

    /// Equality operator
    bool operator==(const NiftiImage &other) const;

    /// Equality operator
    bool operator!=(const NiftiImage &other) const;

    /// Addition operator
    NiftiImage operator+(const NiftiImage&) const;

    /// Subtraction operator
    NiftiImage operator-(const NiftiImage&) const;

    /// Addition operator
    NiftiImage operator+(const float&) const;

    /// Subtraction operator
    NiftiImage operator-(const float&) const;

    /// Multiply image
    NiftiImage operator*(const float &value) const;

    /// Is the image initialised? (Should unless default constructor was used.)
    bool is_initialised() const { return (_nifti_image ? true : false); }

    /// Get image as nifti
    std::shared_ptr<nifti_image> get_raw_nifti_sptr() const;

    /// Save to file
    void save_to_file(const std::string &filename) const;

    /// Get max
    float get_max() const;

    /// Get min
    float get_min() const;

    /// Get element
    float get_element(const int idx[]) const;

    /// Get sum
    float get_sum() const;

    /// Fill
    void fill(const float &v);

    /// Get norm
    float get_norm(const NiftiImage&) const;

    /// Deep copy
    NiftiImage deep_copy() const;

    /// Get number of voxels
    void get_dimensions(int dims[8]) const;

    /// Change image datatype
    template<typename newType>
    void change_datatype();

    /// Get datatype
    std::string get_datatype() const;

    /// Dump header info
    void dump_header() const;

    /// Dump multiple header info
    static void dump_headers(const std::vector<sirf::NiftiImage> &ims);

protected:

    enum NiftiImageType { _general, _3D, _3DTensor, _3DDisp, _3DDef};

    /// Image data as a nifti object
    std::shared_ptr<nifti_image>  _nifti_image;

    /// Check dimensions. Don't require anything for this class.
    void check_dimensions(const enum NiftiImageType image_type = _general);
};
}

#endif
