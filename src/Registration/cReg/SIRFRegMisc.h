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

namespace sirf {
class SIRFRegTransformation;
class SIRFRegTransformationDeformation;
class SIRFImageDataDeformation;
}

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
    void save_multicomponent_nifti_image_split_xyz(std::shared_ptr<nifti_image> input_sptr, const std::string &filename);

    /// Copy nifti image
    void copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr);

    /// Flip multicomponent image along a given axis
    void flip_multicomponent_image(sirf::SIRFImageDataDeformation &im, int dim);

    /// Get cpp from transformation matrix
#if NIFTYREG_VER_1_3
    void get_cpp_from_transformation_matrix(std::shared_ptr<nifti_image> &cpp_sptr, const mat44 &TM_sptr, const std::shared_ptr<nifti_image> &warped_sptr);
#endif
    /// Get def from cpp
    void get_def_from_cpp(sirf::SIRFImageDataDeformation &def, const sirf::SIRFImageDataDeformation &cpp, const sirf::SIRFImageData &ref);

    /// Convert from deformation to displacement field image
    void convert_from_def_to_disp(sirf::SIRFImageDataDeformation &im);

    /// Convert from displacement to deformation field image
    void convert_from_disp_to_def(sirf::SIRFImageDataDeformation &im);

    /// Multiply image
    void multiply_image(sirf::SIRFImageData &output, const sirf::SIRFImageData &input, const float &value);

    /// Do nifti image metadatas match?
    bool do_nifti_image_metadata_match(const sirf::SIRFImageData &im1, const sirf::SIRFImageData &im2);

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
    bool do_nifti_images_match(const sirf::SIRFImageData &im1, const sirf::SIRFImageData &im2, const float accuracy_percentage_of_max/* = 0.F*/);

    /// Do arrays match?
    template<typename T>
    bool do_arrays_match(const sirf::SIRFImageData &im1, const sirf::SIRFImageData &im2, const float accuracy_percentage_of_max)
    {
        if(!im1.is_initialised())
            throw std::runtime_error("do_arrays_match: Image 1 not initialised.");
        if(!im2.is_initialised())
            throw std::runtime_error("do_arrays_match: Image 2 not initialised.");

        // Subtract images
        sirf::SIRFImageData sub = im1 - im2;

        // Get absolute of difference
        T *data = static_cast<T*>(sub.get_raw_nifti_sptr()->data);
        for (unsigned i=0; i<sub.get_raw_nifti_sptr()->nvox; ++i)
            data[i] = fabs(data[i]);

        // Get sum of abs difference
        float sum = sub.get_sum();

        // Get average
        float avg = sum / sub.get_raw_nifti_sptr()->nvox;

        // Calculate required accuracy
        float max1 = im1.get_max();
        float max2 = im2.get_max();
        float epsilon = (max1+max2)/2.F;
        epsilon *= accuracy_percentage_of_max / 100.F;

        if (avg > epsilon) {
            std::cout << "\nImages are not equal.\n";
            std::cout << "\tmax1                              = " << max1 << "\n";
            std::cout << "\tmax2                              = " << max2 << "\n";
            std::cout << "\tmin1                              = " << im1.get_min() << "\n";
            std::cout << "\tmin2                              = " << im2.get_min() << "\n";
            std::cout << "\tpercentage required               = " << accuracy_percentage_of_max << "%\n";
            std::cout << "\tepsilon (avg(max1,max2)/accuracy) = " << epsilon << "\n";
            std::cout << "\tavg abs. diff.                    = " << avg << "\n";
            return false;
        }
        return true;
    }

    /// Get array max
    template<typename T>
    float get_array_max(const sirf::SIRFImageData &im)
    {
        if(!im.is_initialised())
            throw std::runtime_error("get_array_max: Image not initialised.");

        // Check sizes
        if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
            throw std::runtime_error("get_array_max: Datatype does not match desired cast type (" + std::to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        // Get data
        T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);
        return *std::max_element(data, data + im.get_raw_nifti_sptr()->nvox);
    }

    /// Get array min
    template<typename T>
    float get_array_min(const sirf::SIRFImageData &im)
    {
        if(!im.is_initialised())
            throw std::runtime_error("get_array_min: Image not initialised.");

        // Check sizes
        if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
            throw std::runtime_error("get_array_min: Datatype does not match desired cast type (" + std::to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        // Get data
        T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);
        return *std::min_element(data, data + im.get_raw_nifti_sptr()->nvox);
    }

    /// Get array sum
    template<typename T>
    float get_array_sum(const sirf::SIRFImageData &im)
    {
        if(!im.is_initialised())
            throw std::runtime_error("get_array_min: Image not initialised.");

        // Check sizes
        if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
            throw std::runtime_error("get_array_min: Datatype does not match desired cast type (" + std::to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        // Get data
        T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);
        T sum = T(0);
        for (int i=0; i<im.get_raw_nifti_sptr()->nvox; ++i)
            sum += data[i];
        return sum;
    }

    /// Get 3D array element
    template<typename T>
    float get_array_element(const sirf::SIRFImageData &im, int x, int y=0, int z=0, int t=0, int u=0, int v=0, int w=0)
    {
        if(!im.is_initialised())
            throw std::runtime_error("get_3D_array_element: Image not initialised.");

        // Check sizes
        if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
            throw std::runtime_error("get_array_element: Datatype does not match desired cast type (" + std::to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        // Check the point is in bounds
        const int &nx = im.get_raw_nifti_sptr()->nx;
        const int &ny = im.get_raw_nifti_sptr()->ny;
        const int &nz = im.get_raw_nifti_sptr()->nz;
        const int &nt = im.get_raw_nifti_sptr()->nt;
        const int &nu = im.get_raw_nifti_sptr()->nu;
        const int &nv = im.get_raw_nifti_sptr()->ny;
        const int &nw = im.get_raw_nifti_sptr()->nw;
        if(x<0 || x>=nx ||
           y<0 || y>=ny ||
           z<0 || z>=nz ||
           t<0 || t>=nt ||
           u<0 || u>=nu ||
           v<0 || v>=nv ||
           w<0 || w>=nw)
            throw std::runtime_error("get_array_element: Element out of bounds");

        // Get data
        T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);

        std::cout << "\nBe careful, I made this quickly for debugging and haven't thought about data order.\n";

        x *= nw*nv*nu*nt*nz*ny;
        y *= nw*nv*nu*nt*nz;
        z *= nw*nv*nu*nt;
        t *= nw*nv*nu;
        u *= nw*nv;
        v *= nw;
        return data[x+y+z+t+u+v+w];
    }

    /// Sum arrays
    template<typename T>
    sirf::SIRFImageData sum_arrays(const sirf::SIRFImageData &im1, const sirf::SIRFImageData &im2)
    {
        if(!im1.is_initialised())
            throw std::runtime_error("sum_arrays: Image 1 not initialised.");
        if(!im2.is_initialised())
            throw std::runtime_error("sum_arrays: Image 2 not initialised.");
        if(!do_nifti_image_metadata_match(im1,im2))
            throw std::runtime_error("sum_arrays: Cannot add images as metadata does not match.");
        // Check sizes
        if (im1.get_raw_nifti_sptr()->nbyper != sizeof(T))
            throw std::runtime_error("sum_arrays: Datatype of image 1 does not match desired cast type (" + std::to_string(im1.get_raw_nifti_sptr()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        sirf::SIRFImageData result = im1.deep_copy();

        // Get data
        T *im2_data = static_cast<T*>(im2.get_raw_nifti_sptr()->data);
        T *res_data = static_cast<T*>(result.get_raw_nifti_sptr()->data);

        for (unsigned i=0; i<im1.get_raw_nifti_sptr()->nvox; ++i) {
            res_data[i] += im2_data[i];
        }
        return result;
    }

    /// Subtract arrays
    template<typename T>
    sirf::SIRFImageData sub_arrays(const sirf::SIRFImageData &im1, const sirf::SIRFImageData &im2)
    {
        if(!im1.is_initialised())
            throw std::runtime_error("sub_arrays: Image 1 not initialised.");
        if(!im2.is_initialised())
            throw std::runtime_error("sub_arrays: Image 2 not initialised.");
        if(!do_nifti_image_metadata_match(im1,im2))
            throw std::runtime_error("sub_arrays: Cannot add images as metadata does not match.");
        // Check sizes
        if (im1.get_raw_nifti_sptr()->nbyper != sizeof(T))
            throw std::runtime_error("sub_arrays: Datatype of image 1 does not match desired cast type (" + std::to_string(im1.get_raw_nifti_sptr()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        sirf::SIRFImageData result = im1.deep_copy();

        // Get data
        T *im2_data = static_cast<T*>(im2.get_raw_nifti_sptr()->data);
        T *res_data = static_cast<T*>(result.get_raw_nifti_sptr()->data);

        for (unsigned i=0; i<im1.get_raw_nifti_sptr()->nvox; ++i) {
            res_data[i] -= im2_data[i];
        }
        return result;
    }

    /// Dump info of nifti image
    void dump_nifti_info(const std::string &im_filename);

    /// Dump info of nifti image
    void dump_nifti_info(const sirf::SIRFImageData &im);

    /// Dump info of multiple nifti images
    void dump_nifti_info(const std::vector<sirf::SIRFImageData> &ims);

    /// Dump nifti element
    template<typename T>
    void dump_nifti_element(const std::vector<sirf::SIRFImageData> &ims, const std::string &name, const T &call_back)
    {
        std::string header = name + ": ";
        std::cout << "\t" << std::left << std::setw(19) << header;
        for(int i=0; i<ims.size(); i++)
            std::cout << std::setw(19) << ims[i].get_raw_nifti_sptr().get()->*call_back;
        std::cout << "\n";
    }

    /// Dump nifti element
    template<typename T>
    void dump_nifti_element(const std::vector<sirf::SIRFImageData> &ims, const std::string &name, const T &call_back, const unsigned num_elems)
    {
        for(int i=0; i<num_elems; i++) {
            std::string header = name + "[" + std::to_string(i) + "]: ";
            std::cout << "\t" << std::left << std::setw(19) << header;
            for(unsigned j=0; j<ims.size(); j++)
                std::cout << std::setw(19) << (ims[j].get_raw_nifti_sptr().get()->*call_back)[i];
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
    void change_datatype(const sirf::SIRFImageData &image)
    {
        reg_tools_changeDatatype<T>(image.get_raw_nifti_sptr().get());
    }

    /// Fill array with single value.
    template<typename T>
    void fill_array(const sirf::SIRFImageData &im, const float &v)
    {
        if(!im.is_initialised())
            throw std::runtime_error("fill_array: Image not initialised.");

        // Check sizes
        if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
            throw std::runtime_error("fill_array: Datatype does not match desired cast type (" + std::to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + std::to_string(sizeof(T)) + ").");

        // Get data
        T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);
        for (unsigned i=0; i<im.get_raw_nifti_sptr()->nvox; i++) data[i] = T(v);
    }

    /// Compose multiple transformations into single deformation field
    void compose_transformations_into_single_deformation(sirf::SIRFRegTransformationDeformation &def, const std::vector<sirf::SIRFRegTransformation*> &transformations, const sirf::SIRFImageData &ref);

    /// Compose multiple transformations into single deformation field
    void compose_transformations_into_single_deformation(sirf::SIRFRegTransformationDeformation &def, const std::vector<std::shared_ptr<sirf::SIRFRegTransformation> > &transformations, const sirf::SIRFImageData &ref);

    /// Get identity matrix
    mat44 get_identity_matrix();
}

#endif
