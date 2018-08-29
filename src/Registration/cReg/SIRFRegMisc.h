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

#include <iostream>
#include <iomanip>
#include <nifti1_io.h>
#include <boost/filesystem.hpp>
#include <_reg_tools.h>
#include "SIRFImageData.h"

class SIRFImageDataDeformation;

namespace SIRFRegMisc {

    /// Open nifti image
    void open_nifti_image(std::shared_ptr<nifti_image> &image, const boost::filesystem::path filename);

    /// Save nifti image
    void save_nifti_image(nifti_image *image, const std::string filename);

    /// Save nifti image
    void save_nifti_image(std::shared_ptr<nifti_image> image, const std::string filename);

    //! Split multi-component nifti image
    /*! Assume that the multi-component aspect is in dim[5].
        Split the image into that many single-component images
    */
    std::vector<std::shared_ptr<nifti_image> >
        split_multicomponent_nifti_image(std::shared_ptr<nifti_image> input_sptr);

    /// Save a multicomponent nifti image
    void save_multicomponent_nifti_image(std::shared_ptr<nifti_image> input_sptr, const std::string &filename, const bool &split_xyz);

    /// Copy nifti image
    void copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr);

    /// Flip multicomponent image along a given axis
    void flip_multicomponent_image(SIRFImageDataDeformation &im, int dim);

    /// Get cpp from transformation matrix
#if NIFTYREG_VER_1_3
    void get_cpp_from_transformation_matrix(std::shared_ptr<nifti_image> &cpp_sptr, const mat44 &TM_sptr, const std::shared_ptr<nifti_image> &warped_sptr);
#endif
    /// Get def from cpp
    void get_def_from_cpp(SIRFImageDataDeformation &def, const SIRFImageDataDeformation &cpp, const SIRFImageData &ref);

    /// Convert from deformation to displacement field image
    void convert_from_def_to_disp(SIRFImageDataDeformation &im);

    /// Convert from displacement to deformation field image
    void convert_from_disp_to_def(SIRFImageDataDeformation &im);

    /// Multiply image
    void multiply_image(SIRFImageData &output, const SIRFImageData &input, const float &value);

    /// Do nifti image metadatas match?
    bool do_nifti_image_metadata_match(const SIRFImageData &im1, const SIRFImageData &im2);

    /// Do nifti image metadata elements match?
    template<typename T>
    bool do_nifti_image_metadata_elements_match(const std::string &name, const T &elem1, const T &elem2)
    {
        if(fabs(elem1-elem2) < 1.e-7F)
            return true;
        std::cout << "mismatch in " << name << " , (values: " <<  elem1 << " and " << elem2 << ")\n";
        return false;
    }

    /// Do nifti image metadata elements match?
    bool do_nifti_image_metadata_elements_match(const std::string &name, const mat44 &elem1, const mat44 &elem2);

    /// Do nifti images match?
    bool do_nifti_images_match(const SIRFImageData &im1, const SIRFImageData &im2, const float accuracy_percentage_of_max/* = 0.F*/);

    /// Do arrays match?
    template<typename T>
    bool do_arrays_match(const SIRFImageData &im1, const SIRFImageData &im2, const float accuracy_percentage_of_max)
    {
        if(!im1.is_initialised())
            throw std::runtime_error("do_arrays_match: Image 1 not initialised.");
        if(!im2.is_initialised())
            throw std::runtime_error("do_arrays_match: Image 2 not initialised.");

        // Check sizes
        if (im1.get_image_as_nifti()->nbyper != im2.get_image_as_nifti()->nbyper)
            throw std::runtime_error("do_arrays_match: Images are not of same datatype.");
        if (im1.get_image_as_nifti()->nbyper != sizeof(T))
            throw std::runtime_error("do_arrays_match: Datatype does not match desired cast type (" + std::to_string(im1.get_image_as_nifti()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        // Get data
        T *data1 = static_cast<T*>(im1.get_image_as_nifti()->data);
        T *data2 = static_cast<T*>(im2.get_image_as_nifti()->data);

        // Calculate required accuracy
        float max1 = im1.get_max();
        float max2 = im2.get_max();
        float epsilon = max1 > max2 ? max1 : max2;
        epsilon *= accuracy_percentage_of_max;

        for (unsigned i=0; i<im1.get_image_as_nifti()->nvox; ++i)
            if (fabs(data1[i]-data2[i]) > epsilon) {
                std::cout << "\nMismatch in index " << i << " (" << data1[i] << " versus " << data2[i] << ").\n";
                return false;
            }
        return true;
    }

    /// Get array max
    template<typename T>
    float get_array_max(const SIRFImageData &im)
    {
        if(!im.is_initialised())
            throw std::runtime_error("get_array_max: Image not initialised.");

        // Check sizes
        if (im.get_image_as_nifti()->nbyper != sizeof(T))
            throw std::runtime_error("get_array_max: Datatype does not match desired cast type (" + std::to_string(im.get_image_as_nifti()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        // Get data
        T *data = static_cast<T*>(im.get_image_as_nifti()->data);
        return *std::max_element(data, data + im.get_image_as_nifti()->nvox);
    }

    /// Get array min
    template<typename T>
    float get_array_min(const SIRFImageData &im)
    {
        if(!im.is_initialised())
            throw std::runtime_error("get_array_min: Image not initialised.");

        // Check sizes
        if (im.get_image_as_nifti()->nbyper != sizeof(T))
            throw std::runtime_error("get_array_min: Datatype does not match desired cast type (" + std::to_string(im.get_image_as_nifti()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        // Get data
        T *data = static_cast<T*>(im.get_image_as_nifti()->data);
        return *std::min_element(data, data + im.get_image_as_nifti()->nvox);
    }

    /// Get 3D array element
    template<typename T>
    float get_3D_array_element(const SIRFImageData &im, const int x, const int y, const int z)
    {
        if(!im.is_initialised())
            throw std::runtime_error("get_3D_array_element: Image not initialised.");

        // Check sizes
        if (im.get_image_as_nifti()->nbyper != sizeof(T))
            throw std::runtime_error("get_3D_array_element: Datatype does not match desired cast type (" + std::to_string(im.get_image_as_nifti()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        int nx = im.get_image_as_nifti()->nx;
        int ny = im.get_image_as_nifti()->nz;
        int nz = im.get_image_as_nifti()->ny;
        if(x<0 || x>=nx || y<0 || y>=ny || z<0 || z>=nz)
            throw std::runtime_error("get_3D_array_element: Element out of bounds");

        // Get data
        T *data = static_cast<T*>(im.get_image_as_nifti()->data);

        std::cout << "\nBe careful, I made this quickly for debugging and haven't thought about data order."
                     " You might have to switch x and z.\n";

        return data[x*ny*nz + y*nz + z];
    }

    /// Dump info of nifti image
    void dump_nifti_info(const std::string &im_filename);

    /// Dump info of nifti image
    void dump_nifti_info(const SIRFImageData &im);

    /// Dump info of multiple nifti images
    void dump_nifti_info(const std::vector<SIRFImageData> &ims);

    /// Dump nifti element
    template<typename T>
    void dump_nifti_element(const std::vector<SIRFImageData> &ims, const std::string &name, const T &call_back)
    {
        std::cout << "\t" << std::left << std::setw(19) << name << ": ";
        for(int i=0; i<ims.size(); i++)
            std::cout << std::setw(19) << ims[i].get_image_as_nifti().get()->*call_back;
        std::cout << "\n";
    }

    /// Dump nifti element
    template<typename T>
    void dump_nifti_element(const std::vector<SIRFImageData> &ims, const std::string &name, const T &call_back, const unsigned num_elems)
    {
        for(int i=0; i<num_elems; i++) {
            std::cout << "\t" << name << "[" << i << "]:\t\t   ";
            for(unsigned j=0; j<ims.size(); j++)
                std::cout << std::setw(19) << (ims[j].get_image_as_nifti().get()->*call_back)[i];
            std::cout << "\n";
        }
    }

    /// Save transformation matrix to file
    void save_transformation_matrix(const mat44 &transformation_matrix, const std::string &filename);

    /// Read transformation matrix from file
    void open_transformation_matrix(mat44 &transformation_matrix, const std::string &filename);

    /// Print mat44
    void print_mat44(const mat44 &mat);

    /// Print multiple mat44
    void print_mat44(const std::vector<mat44> &mats);

    /// Do mat44 match?
    bool do_mat44_match( const mat44& mat1, const mat44& mat2 );

    /// Mat44 multiplier
    mat44 multiply_mat44(const mat44 &x, const mat44 &y);

    /// Convert type (performs deep copy)
    template<typename T>
    void change_datatype(const SIRFImageData &image)
    {
        reg_tools_changeDatatype<T>(image.get_image_as_nifti().get());
    }
}

#endif
