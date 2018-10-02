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
\brief Generic tools (e.g., opening files)

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFMISC_H_
#define _SIRFMISC_H_

#include <boost/filesystem.hpp>
#include <nifti1_io.h>
#include <sstream>
#include "NiftiImage3D.h"

namespace sirf {
class NiftiImage3DTensor;
class NiftiImage3DDisplacement;
class NiftiImage3DDeformation;
class SIRFRegTransformation;
class SIRFRegMat44;
}

namespace SIRFRegMisc {

    /// Open nifti image
    void open_nifti_image(std::shared_ptr<nifti_image> &image, const boost::filesystem::path &filename);

    /// Save nifti image
    void save_nifti_image(const sirf::NiftiImage &image, const std::string &filename);

    //! Split multi-component nifti image
    /*! Assume that the multi-component aspect is in dim[5].
        Split the image into that many single-component images
    */
    std::vector<sirf::NiftiImage3D> split_multicomponent_nifti_image(const sirf::NiftiImage3DTensor &input);

    /// Save a multicomponent nifti image as x-, y-, z-components
    void save_multicomponent_nifti_image_split_xyz(const sirf::NiftiImage3DTensor &input, const std::string &filename);

    /// Save a multicomponent nifti image as x-, y-, z-components
    void save_multicomponent_nifti_image_split_xyz(const sirf::NiftiImage3DTensor &input, const std::string &filename_x, const std::string &filename_y, const std::string &filename_z);

    /// Copy nifti image
    void copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr);

    /// Flip multicomponent image along a given axis
    void flip_multicomponent_image(sirf::NiftiImage3DTensor &im, int dim);

    /// Get cpp from transformation matrix
#if NIFTYREG_VER_1_3
    void get_cpp_from_transformation_matrix(std::shared_ptr<nifti_image> &cpp_sptr, const mat44 &TM_sptr, const std::shared_ptr<nifti_image> &warped_sptr);
#endif
    /// Get def from cpp
    void get_def_from_cpp(sirf::NiftiImage3DDeformation &def, const sirf::NiftiImage3DTensor &cpp, const sirf::NiftiImage3D &ref);

    /// Convert from deformation to displacement field image
    void convert_from_def_to_disp(sirf::NiftiImage3DDisplacement &disp, const sirf::NiftiImage3DDeformation &def);

    /// Convert from displacement to deformation field image
    void convert_from_disp_to_def(sirf::NiftiImage3DDeformation &def, const sirf::NiftiImage3DDisplacement &disp);

    /// Do nifti image metadatas match?
    bool do_nifti_image_metadata_match(const sirf::NiftiImage &im1, const sirf::NiftiImage &im2);

    /// Do nifti image metadata elements match?
    template<typename T>
    bool do_nifti_image_metadata_elements_match(const std::string &name, const T &elem1, const T &elem2);

    /// Do nifti image metadata elements match?
    bool do_nifti_image_metadata_elements_match(const std::string &name, const mat44 &elem1, const mat44 &elem2);

    /// Do nifti images match?
    bool do_nifti_images_match(const sirf::NiftiImage &im1, const sirf::NiftiImage &im2, const float accuracy_percentage_of_max/* = 0.F*/);

    /// Do arrays match?
    template<typename T>
    bool do_arrays_match(const sirf::NiftiImage &im1, const sirf::NiftiImage &im2, const float &required_accuracy_compared_to_max);

    /// Get array max
    template<typename T>
    float get_array_max(const sirf::NiftiImage &im);
    /// Get array min
    template<typename T>
    float get_array_min(const sirf::NiftiImage &im);

    /// Get array sum
    template<typename T>
    float get_array_sum(const sirf::NiftiImage &im);

    /// Get 3D array element. idx can have up to 7 dims
    template<typename T>
    float get_array_element(const sirf::NiftiImage &im, const int idx[7]);

    /// Sum arrays
    template<typename T>
    sirf::NiftiImage sum_arrays(const sirf::NiftiImage &im1, const sirf::NiftiImage &im2);

    /// Subtract arrays
    template<typename T>
    sirf::NiftiImage sub_arrays(const sirf::NiftiImage &im1, const sirf::NiftiImage &im2);

    /// Array norm
    template<typename T>
    float arrays_norm(const sirf::NiftiImage &im1, const sirf::NiftiImage &im2);

    /// Dump info of multiple nifti images
    void dump_headers_actual(const std::vector<sirf::NiftiImage> &ims);

    /// Dump nifti element
    template<typename T>
    void dump_nifti_element(const std::vector<sirf::NiftiImage> &ims, const std::string &name, const T &call_back);

    /// Dump nifti element
    template<typename T>
    void dump_nifti_element(const std::vector<sirf::NiftiImage> &ims, const std::string &name, const T &call_back, const unsigned num_elems);

    /// Convert type (performs deep copy)
    template<typename newType, typename oldType>
    void change_datatype1(const sirf::NiftiImage &image)
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
            ss << "SIRFRegMisc::change_datatype not implemented for your new data type: ";
            ss << typeid (newType).name();
            ss << " (bytes per voxel: ";
            ss << sizeof(newType) << ").";
            throw std::runtime_error(ss.str());
        }

        im->nbyper = sizeof(newType);
        im->data = static_cast<void*>(calloc(im->nvox,sizeof(newType)));
        newType *dataPtr = static_cast<newType*>(im->data);
        for (size_t i = 0; i < im->nvox; i++)
           dataPtr[i] = static_cast<newType>(originalArray[i]);

        free(originalArray);
        return;
    }

    /// Fill array with single value.
    template<typename T>
    void fill_array(const sirf::NiftiImage &im, const float &v);

    /// Crop image
    template<typename T>
    void crop_image(sirf::NiftiImage &image, const int min_index[7], const int max_index[7]);
}

#endif
