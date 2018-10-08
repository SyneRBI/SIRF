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

    /// Copy constructor
    NiftiImage(const NiftiImage& to_copy);

    /// Assignment
    NiftiImage& operator=(const NiftiImage& to_copy);

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
    NiftiImage operator*(const float&) const;

    /// Access data element via 1D index (const)
    float operator()(const int index) const;

    /// Access data element via 1D index
    float &operator()(const int index);

    /// Access data element via 7D index (const)
    float operator()(const int index[7]) const;

    /// Access data element via 7D index
    float &operator()(const int index[7]);

    /// Is the image initialised? (Should unless default constructor was used.)
    bool is_initialised() const { return (_nifti_image && _data && _nifti_image->datatype == DT_FLOAT32 ? true : false); }

    /// Get image as nifti as const
    std::shared_ptr<const nifti_image> get_raw_nifti_sptr() const;

    /// Get image as nifti
    std::shared_ptr<nifti_image> get_raw_nifti_sptr();

    /// Save to file. Templated so the user can choose the datatype they save to. This defaults
    /// to -1, which is the original datatype of that image (stored as _original_datatype).
    void save_to_file(const std::string &filename, const int datatype = -1) const;

    /// Get max
    float get_max() const;

    /// Get min
    float get_min() const;

    /// Get mean
    float get_mean() const;

    /// Get element
    float get_element(const int idx[7]) const;

    /// Get sum
    float get_sum() const;

    /// Fill
    void fill(const float &v);

    /// Get norm
    float get_norm(const NiftiImage&) const;

    /// Deep copy
    NiftiImage deep_copy() const;

    /// Get number of voxels
    const int* get_dimensions() const;

    /// Print header info
    void print_header() const;

    /// Print multiple header info
    static void print_headers(const std::vector<sirf::NiftiImage> &ims);

    /// Crop
    void crop(const int min_index[7], const int max_index[7]);

    /// get 1D index from ND index
    int get_1D_index(const int idx[7]) const;

    /// Get original datatype
    int get_original_datatype() const { return _original_datatype; }

    /// Check if the norms of two images are equal to a given accuracy.
    bool are_equal_to_given_accuracy(const NiftiImage &im2, const float required_accuracy_compared_to_max) const;

    /// Point is in bounds?
    bool is_in_bounds(const int index[7]) const;

    /// Point is in bounds?
    bool is_in_bounds(const int index) const;

    /// Images are same size
    bool is_same_size(const NiftiImage &im) const;

protected:

    enum NiftiImageType { _general, _3D, _3DTensor, _3DDisp, _3DDef};

    enum MathsType { add, sub, mul };

    /// Image data as a nifti object
    std::shared_ptr<nifti_image>  _nifti_image;

    /// Data
    float *_data = NULL;

    /// Original datatype
    int _original_datatype = -1;

    /// Check dimensions. Don't require anything for this class.
    void check_dimensions(const enum NiftiImageType image_type = _general);

    /// Set up datatype. Set to float if not already, store the original type.
    void set_up_data(const int original_datatype);

    /// Add, subract image from another
    NiftiImage maths(const NiftiImage& c, const MathsType type) const;

    /// Add, subract, multiply value to image
    NiftiImage maths(const float val, const MathsType type) const;

private:

    /// Change image datatype with int
    void change_datatype(const int datatype);
};
}

#endif
