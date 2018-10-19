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
#include "NiftiImageData3D.h"

namespace sirf {
class NiftiImageData3DTensor;
class NiftiImageData3DDisplacement;
class NiftiImageData3DDeformation;
class SIRFRegTransformation;
class SIRFRegAffineTransformation;
}

namespace SIRFRegMisc {

    /// Open nifti image
    void open_nifti_image(std::shared_ptr<nifti_image> &image, const boost::filesystem::path &filename);

    /// Save nifti image. image is not const because filename gets set during the process.
    void save_nifti_image(sirf::NiftiImageData &image, const std::string &filename);

    /// Copy nifti image
    void copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr);

    /// Do nifti image metadatas match?
    bool do_nifti_image_metadata_match(const sirf::NiftiImageData &im1, const sirf::NiftiImageData &im2);

    /// Do nifti image metadata elements match?
    template<typename T>
    bool do_nifti_image_metadata_elements_match(const std::string &name, const T &elem1, const T &elem2);

    /// Do nifti image metadata elements match?
    bool do_nifti_image_metadata_elements_match(const std::string &name, const mat44 &elem1, const mat44 &elem2);

    /// Dump info of multiple nifti images
    void dump_headers(const std::vector<sirf::NiftiImageData> &ims);

    /// Dump nifti element
    template<typename T>
    void dump_nifti_element(const std::vector<sirf::NiftiImageData> &ims, const std::string &name, const T &call_back);

    /// Dump nifti element
    template<typename T>
    void dump_nifti_element(const std::vector<sirf::NiftiImageData> &ims, const std::string &name, const T &call_back, const unsigned num_elems);

    /// Change datatype. Templated for desired type. Figures out what current type is then calls doubley templated function below.
    template<typename newType>
    void change_datatype(sirf::NiftiImageData &im);

    /// Convert type (performs deep copy)
    template<typename newType, typename oldType>
    void change_datatype(sirf::NiftiImageData &image)
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
}

#endif
