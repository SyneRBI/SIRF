/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2019 University College London

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

#pragma once

#include <nifti1_io.h>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <sstream>
#include "sirf/common/ANumRef.h"
#include "sirf/common/ImageData.h"
#include <_reg_tools.h>

namespace sirf {

/*!
\file
\ingroup Registration
\brief Base class for SIRF's nifti image data.

This is a wrapper around the basic nifti_image (stored as a sptr),
with extra functionality.

This is the general form, and any number of dimensions are allowed. 
This is contrary to the derived classes, which have specific requirements
(e.g., NiftiImageData3DDeformation requires 5 dimensions, of which the 4th (time) == 1).

To move between different images types (e.g., STIRImageData and MRImageData), 
we need to know the patient's position relative to the scanner. For this, only
qform_code == 1 will suffice.

qform/sform
The qform code describes "scanner-anatomical" coordinates, whereas the
sform code describes the location of the voxels in some standard space.

For qform > 0, the origin of coordinates would generally be whatever
the scanner origin is; for example, in MRI, (0,0,0) is the center
of the gradient coil.

For sform > 0, the origin of coordinates would depend on the value
of sform_code; for example, for the Talairach coordinate system,
(0,0,0) corresponds to the Anterior Commissure.

\author Richard Brown
\author CCP PETMR
*/

template<class dataType>
class NiftiImageData : public ImageData
{
public:

    typedef ImageData::Iterator BaseIter;
    typedef ImageData::Iterator_const BaseIter_const;
    class Iterator : public BaseIter {
    public:
        Iterator(dataType *iter) : _iter(iter)
        {}
        Iterator& operator=(const Iterator& iter)
        {
            _iter = iter._iter;
            _ref.copy(iter._ref);
            return *this;
        }
        virtual Iterator& operator++()
        {
            ++_iter;
            return *this;
        }
        virtual bool operator==(const BaseIter& an_iter) const
        {
            const Iterator& iter = (const Iterator&)an_iter;
            return _iter == iter._iter;
        }
        virtual bool operator!=(const BaseIter& an_iter) const
        {
            const Iterator& iter = (const Iterator&)an_iter;
            return _iter != iter._iter;
        }
        virtual FloatRef& operator*()
        {
            dataType& v = *_iter;
            _ref.set_ptr(&v);
            return _ref;
        }
    private:
        dataType *_iter;
        FloatRef _ref;
    };
    class Iterator_const : public BaseIter_const {
    public:
        Iterator_const(const dataType *iter) : _iter(iter)
        {}
        Iterator_const& operator=(const Iterator_const& iter)
        {
            _iter = iter._iter;
            _ref.copy(iter._ref);
            return *this;
        }
        virtual Iterator_const& operator++()
        {
            ++_iter;
            return *this;
        }
        virtual bool operator==(const BaseIter_const& an_iter) const
        {
            const Iterator_const& iter = (const Iterator_const&)an_iter;
            return _iter == iter._iter;
        }
        virtual bool operator!=(const BaseIter_const& an_iter) const
        {
            const Iterator_const& iter = (const Iterator_const&)an_iter;
            return _iter != iter._iter;
        }
        virtual const FloatRef& operator*() const
        {
            const dataType& v = *_iter;
            _ref.set_ptr((void*)&v);
            return _ref;
        }
    private:
        const dataType *_iter;
        mutable FloatRef _ref;
    };

    /// Constructor
    NiftiImageData() {}

    /// Destructor
    virtual ~NiftiImageData() {}

    /// Copy constructor
    NiftiImageData(const NiftiImageData& to_copy);

    /// Assignment
    NiftiImageData& operator=(const NiftiImageData& to_copy);

    /// Copy constructor
    NiftiImageData(const ImageData& to_copy);

    /// Assignment
    NiftiImageData& operator=(const ImageData& to_copy);

    /// Filename constructor
    NiftiImageData(const std::string &filename);

    /// Nifti constructor
    NiftiImageData(const nifti_image &image_nifti);

    /// Construct from array
    template<class inputType>
    NiftiImageData(const inputType * const data, const VoxelisedGeometricalInfo3D &geom, const bool is_tensor = false)
    {
        this->_nifti_image = create_from_geom_info(geom, is_tensor);

        // Set the datatype
        if      (typeid(inputType) == typeid(bool))               this->set_up_data(DT_BINARY);
        else if (typeid(inputType) == typeid(signed char))        this->set_up_data(DT_INT8);
        else if (typeid(inputType) == typeid(signed short))       this->set_up_data(DT_INT16);
        else if (typeid(inputType) == typeid(signed int))         this->set_up_data(DT_INT32);
        else if (typeid(inputType) == typeid(float))              this->set_up_data(DT_FLOAT32);
        else if (typeid(inputType) == typeid(double))             this->set_up_data(DT_FLOAT64);
        else if (typeid(inputType) == typeid(unsigned char))      this->set_up_data(DT_UINT8);
        else if (typeid(inputType) == typeid(unsigned short))     this->set_up_data(DT_UINT16);
        else if (typeid(inputType) == typeid(unsigned int))       this->set_up_data(DT_UINT32);
        else if (typeid(inputType) == typeid(signed long long))   this->set_up_data(DT_INT64);
        else if (typeid(inputType) == typeid(unsigned long long)) this->set_up_data(DT_UINT64);
        else if (typeid(inputType) == typeid(long double))        this->set_up_data(DT_FLOAT128);
        else {
            std::stringstream ss;
            ss << "NiftiImageData constructor from raw array: ";
            ss << typeid (inputType).name();
            ss << " (bytes per voxel: ";
            ss << sizeof(inputType) << ").";
            throw std::runtime_error(ss.str());
        }

        for (unsigned i=0; i<_nifti_image->nvox; ++i)
            this->_data[i] = dataType(data[i]);
    }

    /// Create NiftiImageData from geometrical info
    static std::shared_ptr<nifti_image> create_from_geom_info(const VoxelisedGeometricalInfo3D &geom, const bool is_tensor=false);

    /// Construct NiftiImageData from the real component of a complex SIRF ImageData
    static void construct_NiftiImageData_from_complex_im_real_component(std::shared_ptr<NiftiImageData> &out_sptr, const std::shared_ptr<const ImageData> in_sptr);

    /// Construct NiftiImageData from the imaginary component of a complex SIRF ImageData
    static void construct_NiftiImageData_from_complex_im_imag_component(std::shared_ptr<NiftiImageData> &out_sptr, const std::shared_ptr<const ImageData> in_sptr);

    /// Construct two NiftiImageData from a complex SIRF ImageData
    static void construct_NiftiImageData_from_complex_im(std::shared_ptr<NiftiImageData> &out_real_sptr, std::shared_ptr<NiftiImageData> &out_imag_sptr, const std::shared_ptr<const ImageData> in_sptr);

    /// Equality operator
    bool operator==(const NiftiImageData &other) const;

    /// Equality operator
    bool operator!=(const NiftiImageData &other) const;

    /// Addition operator
    NiftiImageData& operator+=(const NiftiImageData &rhs);

    /// Addition operator
    friend NiftiImageData operator+(NiftiImageData lhs, const NiftiImageData& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    /// Subtraction operator
    NiftiImageData& operator-=(const NiftiImageData &rhs);

    /// Subtraction operator
    friend NiftiImageData operator-(NiftiImageData lhs, const NiftiImageData& rhs)
    {
        lhs -= rhs;
        return lhs;
    }

    /// Addition operator
    NiftiImageData& operator+=(const float);

    /// Addition operator
    friend NiftiImageData operator+(NiftiImageData lhs, const float val)
    {
        lhs += val;
        return lhs;
    }

    /// Subtraction operator
    NiftiImageData& operator-=(const float);

    /// Subtraction operator
    friend NiftiImageData operator-(NiftiImageData lhs, const float val)
    {
        lhs -= val;
        return lhs;
    }

    /// Multiplication operator
    NiftiImageData& operator*=(const float);

    /// Multiplication operator
    friend NiftiImageData operator*(NiftiImageData lhs, const float val)
    {
        lhs *= val;
        return lhs;
    }

    /// Division operator
    NiftiImageData& operator/=(const float);

    /// Division operator
    friend NiftiImageData operator/(NiftiImageData lhs, const float val)
    {
        lhs /= val;
        return lhs;
    }

    /// Access data element via 1D index (const)
    float operator()(const int index) const;

    /// Access data element via 1D index
    float &operator()(const int index);

    /// Access data element via 7D index (const)
    float operator()(const int index[7]) const;

    /// Access data element via 7D index
    float &operator()(const int index[7]);

    /// Is the image initialised? (Should be unless default constructor was used.)
    bool is_initialised() const { return (_nifti_image && _data && _nifti_image->datatype == DT_FLOAT32 ? true : false); }

    /// Get image as nifti as const
    std::shared_ptr<const nifti_image> get_raw_nifti_sptr() const;

    /// Get image as nifti
    std::shared_ptr<nifti_image> get_raw_nifti_sptr();

    /// Save to file. Templated so the user can choose the datatype they save to. This defaults
    /// to -1, which is the original datatype of that image (stored as _original_datatype).
    virtual void write(const std::string &filename, const int datatype) const;

    /// Write
    virtual void write(const std::string &filename) const { this->write(filename,-1); }

    /// Get max
    float get_max() const;

    /// Get min
    float get_min() const;

    /// Get mean
    float get_mean() const;

    /// Get variance
    float get_variance() const;

    /// Get standard deviation
    float get_standard_deviation() const;

    /// Get element
    float get_element(const int idx[7]) const;

    /// Get sum
    float get_sum() const;

    /// Get nan count
    unsigned get_nan_count() const;

    /// Fill
    void fill(const float v);

    /// Get norm
    float get_norm(const NiftiImageData&) const;

    /// Get data size in each dimension
    const int* get_dimensions() const;

    /// Get total number of voxels
    size_t get_num_voxels() const;

    /// Print header info
    void print_header() const;

    /// Print multiple header info
    static void print_headers(const std::vector<const NiftiImageData*> &ims);

    /// Crop. Set to -1 to leave unchanged
    void crop(const int min_index[7], const int max_index[7]);

    /// Pad image with value. Give number of voxels to increase in min and max directions. Set values to -1 to leave unchanged
    void pad(const int min_index[7], const int max_index[7], const dataType val = 0);

    /// get 1D index from ND index
    int get_1D_index(const int idx[7]) const;

    /// Get original datatype
    int get_original_datatype() const { return _original_datatype; }

    /// Check if the norms of two images are equal to a given accuracy.
    static bool are_equal_to_given_accuracy(const std::shared_ptr<const NiftiImageData> &im1_sptr, const std::shared_ptr<const NiftiImageData> &im2_sptr, const float required_accuracy_compared_to_max);

    /// Check if the norms of two images are equal to a given accuracy.
    static bool are_equal_to_given_accuracy(const NiftiImageData &im1, const NiftiImageData &im2, const float required_accuracy_compared_to_max);

    /// Point is in bounds?
    bool is_in_bounds(const int index[7]) const;

    /// Point is in bounds?
    bool is_in_bounds(const int index) const;

    /// Images are same size
    bool is_same_size(const NiftiImageData &im) const;

    /// Do nifti image metadatas match?
    static bool do_nifti_image_metadata_match(const NiftiImageData &im1, const NiftiImageData &im2, bool verbose);

    /// Dump info of multiple nifti images
    static void dump_headers(const std::vector<const NiftiImageData*> &ims);

    /// Dump nifti element
    template<typename T>
    static void dump_nifti_element(const std::vector<const NiftiImageData*> &ims, const std::string &name, const T &call_back);

    /// Dump nifti element
    template<typename T>
    static void dump_nifti_element(const std::vector<const NiftiImageData*> &ims, const std::string &name, const T &call_back, const unsigned num_elems);

    /// Set the voxel spacing. Requires resampling image, and so interpolation order is required.
    /// As per NiftyReg, interpolation_order can be either 0, 1 or 3 meaning nearest neighbor, linear or cubic spline interpolation.
    void set_voxel_spacing(const float factors[3], const int interpolation_order);

    /// Kernel convolution
    void kernel_convolution(const float sigma, NREG_CONV_KERNEL_TYPE conv_type = GAUSSIAN_KERNEL);

    /// Does the image contain any NaNs?
    bool get_contains_nans() const { return (this->get_nan_count() > 0); }

    /// Flip the image along a given axis (Rotation of 180 degrees about axis)
    void flip_along_axis(const unsigned axis);

    /// Mirror the image along a given axis (This will change handedness of image)
    void mirror_along_axis(const unsigned axis);

    /// Inner product of two images.
    dataType get_inner_product(const NiftiImageData &other) const;

    /// Normalise image between 0 and 1
    void normalise_zero_and_one();

    /// Standardise (subtract mean and divide by standard deviation).
    void standardise();

protected:

    enum NiftiImageDataType { _general, _3D, _3DTensor, _3DDisp, _3DDef};

    enum MathsType { add, sub, mul };

    /// Image data as a nifti object
    std::shared_ptr<nifti_image>  _nifti_image;

    /// Data
    float *_data = NULL;

    /// Original datatype
    int _original_datatype = -1;

    /// Check dimensions. Don't require anything for this class.
    void check_dimensions(const enum NiftiImageDataType image_type = _general);

    /// Set up datatype. Set to float if not already, store the original type.
    void set_up_data(const int original_datatype);

    /// Add, subract image from another
    void maths(const NiftiImageData& c, const MathsType type);

    /// Add, subract, multiply value to image
    void maths(const float val, const MathsType type);

    /// Open nifti image
    static void open_nifti_image(std::shared_ptr<nifti_image> &image, const std::string &filename);

    /// Copy nifti image
    static void copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr);

private:

    /// Change image datatype with int. Values can be found in nifti1.h (e.g., #define DT_BINARY 1)
    void change_datatype(const int datatype);

    /// Change datatype. Templated for desired type. Figures out what current type is then calls doubley templated function below.
    template<typename newType>
    static void change_datatype(NiftiImageData<dataType> &im)
    {
        if (im.get_raw_nifti_sptr()->datatype == DT_BINARY)   return change_datatype<newType,bool>              (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_INT8)     return change_datatype<newType,signed char>       (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_INT16)    return change_datatype<newType,signed short>      (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_INT32)    return change_datatype<newType,signed int>        (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT32)  return change_datatype<newType,float>             (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT64)  return change_datatype<newType,double>            (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_UINT8)    return change_datatype<newType,unsigned char>     (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_UINT16)   return change_datatype<newType,unsigned short>    (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_UINT32)   return change_datatype<newType,unsigned int>      (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_INT64)    return change_datatype<newType,signed long long>  (im);
        if (im.get_raw_nifti_sptr()->datatype == DT_UINT64)   return change_datatype<newType,unsigned long long>(im);
        if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT128) return change_datatype<newType,long double>       (im);

        std::stringstream ss;
        ss << "change_datatype not implemented for your data type: ";
        ss << nifti_datatype_string(im.get_raw_nifti_sptr()->datatype);
        ss << " (bytes per voxel: ";
        ss << im.get_raw_nifti_sptr()->nbyper << ").";
        throw std::runtime_error(ss.str());
    }

    /// Convert type (performs deep copy)
    template<typename newType, typename oldType>
    static void change_datatype(NiftiImageData &image)
    {
        // If the two types are equal, nothing to be done.
        if (typeid (newType) == typeid(oldType))
            return;

        nifti_image *im = image.get_raw_nifti_sptr().get();

        // Copy the original array
        oldType *originalArray = static_cast<oldType*>(malloc(im->nvox*im->nbyper));
        memcpy(originalArray, im->data, im->nvox*im->nbyper);
        // Reset image array
        free(im->data);

        // Set the datatype
        if      (typeid(newType) == typeid(bool))               im->datatype = DT_BINARY;
        else if (typeid(newType) == typeid(signed char))        im->datatype = DT_INT8;
        else if (typeid(newType) == typeid(signed short))       im->datatype = DT_INT16;
        else if (typeid(newType) == typeid(signed int))         im->datatype = DT_INT32;
        else if (typeid(newType) == typeid(float))              im->datatype = DT_FLOAT32;
        else if (typeid(newType) == typeid(double))             im->datatype = DT_FLOAT64;
        else if (typeid(newType) == typeid(unsigned char))      im->datatype = DT_UINT8;
        else if (typeid(newType) == typeid(unsigned short))     im->datatype = DT_UINT16;
        else if (typeid(newType) == typeid(unsigned int))       im->datatype = DT_UINT32;
        else if (typeid(newType) == typeid(signed long long))   im->datatype = DT_INT64;
        else if (typeid(newType) == typeid(unsigned long long)) im->datatype = DT_UINT64;
        else if (typeid(newType) == typeid(long double))        im->datatype = DT_FLOAT128;
        else {
            std::stringstream ss;
            ss << "change_datatype not implemented for your new data type: ";
            ss << typeid (newType).name();
            ss << " (bytes per voxel: ";
            ss << sizeof(newType) << ").";
            throw std::runtime_error(ss.str());
        }

        // Set the nbyper and swap size from the datatype
        nifti_datatype_sizes(im->datatype, &im->nbyper, &im->swapsize);

        // Copy data
        im->data = static_cast<void*>(calloc(im->nvox,sizeof(newType)));
        newType *dataPtr = static_cast<newType*>(im->data);
        for (size_t i = 0; i < im->nvox; i++)
           dataPtr[i] = newType(originalArray[i]);

        free(originalArray);
        return;
    }

    // ------------------------------------------------------------------------------ //
    // Pure virtual methods from ImageData
    // ------------------------------------------------------------------------------ //
public:
    /// Clone and return as unique pointer.
    std::unique_ptr<NiftiImageData> clone() const
    {
	return std::unique_ptr<NiftiImageData>(this->clone_impl());
    }
    virtual Iterator& begin()
    {
        _begin.reset(new Iterator(_data));
        return *_begin;
    }
    virtual Iterator_const& begin() const
    {
        _begin_const.reset(new Iterator_const(_data));
        return *_begin_const;
    }
    virtual Iterator& end()
    {
        _end.reset(new Iterator(_data+_nifti_image->nvox));
        return *_end;
    }
    virtual Iterator_const& end() const
    {
        _end_const.reset(new Iterator_const(_data+_nifti_image->nvox));
        return *_end_const;
    }
protected:
    /// Clone helper function. Don't use.
    virtual NiftiImageData* clone_impl() const
    {
	return new NiftiImageData(*this);
    }
    virtual ObjectHandle<DataContainer>* new_data_container_handle() const
    {
        return new ObjectHandle<DataContainer>
            (std::shared_ptr<DataContainer>(new NiftiImageData));
    }
    unsigned int items() const { return 1; }
    virtual void dot      (const DataContainer& a_x, void* ptr) const;
    virtual void axpby    (const void* ptr_a, const DataContainer& a_x, const void* ptr_b, const DataContainer& a_y);
    virtual float norm() const;
    virtual void multiply (const DataContainer& a_x, const DataContainer& a_y);
    virtual void divide   (const DataContainer& a_x, const DataContainer& a_y);
    virtual Dimensions dimensions() const
    {
        Dimensions dim;
        int *d = _nifti_image->dim;
        dim["x"] = d[1];
        dim["y"] = d[2];
        dim["z"] = d[3];
        dim["t"] = d[4];
        dim["u"] = d[5];
        dim["v"] = d[6];
        dim["w"] = d[7];
        return dim;
    }
    /// Set up the geometrical info. Use qform preferentially over sform.
    virtual void set_up_geom_info();
protected:
    mutable std::shared_ptr<Iterator> _begin;
    mutable std::shared_ptr<Iterator> _end;
    mutable std::shared_ptr<Iterator_const> _begin_const;
    mutable std::shared_ptr<Iterator_const> _end_const;
};
}
