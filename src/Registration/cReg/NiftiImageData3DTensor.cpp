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
\brief Class for tensor SIRF image data.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftiImageData3DTensor.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include <boost/format.hpp>

using namespace sirf;

template<class dataType>
NiftiImageData3DTensor<dataType>::NiftiImageData3DTensor(const NiftiImageData3D<dataType> &x, const NiftiImageData3D<dataType> &y, const NiftiImageData3D<dataType> &z)
{
    // Check everything is intialised
    if (!x.is_initialised() || !y.is_initialised() || !z.is_initialised())
        throw std::runtime_error("NiftiImageData3DTensor: x,y,z->tensor: Can't create from separate 3D components, as some are uninitialised.");

    if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(x,y,true))
        throw std::runtime_error("NiftiImageData3DTensor: x,y,z->tensor: x and y components don't match.");
    if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(x,z,true))
        throw std::runtime_error("NiftiImageData3DTensor: x,y,z->tensor: x and z components don't match.");

    // Create a 4D from one of the components
    this->create_from_3D_image(x);
    std::vector<NiftiImageData3D<dataType> > ims{x, y, z};

    // for nu=3, the tensor data is stored last.
    //So memcpy x into first third, y into second third and z into last third
    size_t mem = x.get_raw_nifti_sptr()->nvox;

    for (size_t i=0; i<3; ++i) {
        size_t index = mem*i;
        float *src  = static_cast<float*>(ims[i].get_raw_nifti_sptr()->data);
        memcpy(this->_data+index, src, mem*size_t(ims[i].get_raw_nifti_sptr()->nbyper));
    }
}

template<class dataType>
void NiftiImageData3DTensor<dataType>::create_from_3D_image(const NiftiImageData<dataType> &image)
{
    if (!image.is_initialised())
        throw std::runtime_error("NiftiImageData3DTensor::create_from_3D_image. Input image not initialised.");

    std::shared_ptr<const nifti_image> image_sptr = image.get_raw_nifti_sptr();

    // Calculate deformation field image
    nifti_image *output_ptr;
    output_ptr = nifti_copy_nim_info(image_sptr.get());
    output_ptr->dim[0]=output_ptr->ndim=5;
    output_ptr->dim[1]=output_ptr->nx=image_sptr->nx;
    output_ptr->dim[2]=output_ptr->ny=image_sptr->ny;
    output_ptr->dim[3]=output_ptr->nz=image_sptr->nz;
    output_ptr->dim[4]=output_ptr->nt=1;output_ptr->pixdim[4]=output_ptr->dt=1.0;
    if(image_sptr->nz>1) output_ptr->dim[5]=output_ptr->nu=3;
    else output_ptr->dim[5]=output_ptr->nu=2;
    output_ptr->pixdim[5]=output_ptr->du=1.0;
    output_ptr->dim[6]=output_ptr->nv=1;output_ptr->pixdim[6]=output_ptr->dv=1.0;
    output_ptr->dim[7]=output_ptr->nw=1;output_ptr->pixdim[7]=output_ptr->dw=1.0;
    output_ptr->nvox=unsigned(output_ptr->nx*output_ptr->ny*output_ptr->nz*output_ptr->nt*output_ptr->nu);
    output_ptr->datatype = DT_FLOAT32;
    output_ptr->nbyper = sizeof(float);
    output_ptr->data = static_cast<void *>(calloc(output_ptr->nvox, unsigned(output_ptr->nbyper)));
    output_ptr->intent_code = NIFTI_INTENT_VECTOR;

    this->_nifti_image = std::shared_ptr<nifti_image>(output_ptr, nifti_image_free);

    this->set_up_data(DT_FLOAT32);
}

template<class dataType>
std::shared_ptr<NiftiImageData<dataType> > NiftiImageData3DTensor<dataType>::get_tensor_component(const int component) const
{
    // Check in range
    if (component < 0 || component >= this->_nifti_image->nu)
        throw std::runtime_error("NiftiImageData3DTensor<dataType>::get_tensor_component: requested component out of range.");

    // Crop image
    std::shared_ptr<NiftiImageData<dataType> > image_sptr = this->clone();
    int index[7] = {-1, -1, -1, -1, component, -1, -1};
    image_sptr->crop(index,index);

    // Intent code is no longer vector
    image_sptr->get_raw_nifti_sptr()->intent_code = NIFTI_INTENT_NONE;
    image_sptr->get_raw_nifti_sptr()->intent_p1 = -1;

    return image_sptr;
}

template<class dataType>
void NiftiImageData3DTensor<dataType>::write_split_xyz_components(const std::string &filename_pattern, const int datatype) const
{
    // Check that the disp image exists
    if (!this->is_initialised())
        throw std::runtime_error("Error, cannot write " + filename_pattern + " to file because it has not been initialised.");

    // Check that filename isn't blank
    if (filename_pattern.empty())
        throw std::runtime_error("Error, cannot write " + filename_pattern + " to file because filename is blank.");

    // Check it contains the %s format
    size_t pos = filename_pattern.find("%s");
    if (pos == filename_pattern.npos)
        throw std::runtime_error("Filename (" + filename_pattern + ") should be given in boost format (e.g., output_%s.nii)");

    std::vector<std::string> filenames(3);
    for (unsigned i=0; i<3; ++i)
        filenames[i] = filename_pattern;

    filenames[0].replace(pos, 2, "x");
    filenames[1].replace(pos, 2, "y");
    filenames[2].replace(pos, 2, "z");

    this->write_split_xyz_components(filenames[0], filenames[1], filenames[2], datatype);
}

template<class dataType>
void NiftiImageData3DTensor<dataType>::write_split_xyz_components(const std::string &filename_x, const std::string &filename_y, const std::string &filename_z, const int datatype) const
{
    int min_index[7], max_index[7];
    for (int i=0; i<7; ++i) {
        min_index[i] = 0;
        max_index[i] = this->_nifti_image->dim[i+1] - 1;
    }

    for (int i=0; i<3; ++i) {

        std::shared_ptr<NiftiImageData<dataType> > image_sptr =
                this->get_tensor_component(i);

        if      (i == 0) image_sptr->write(filename_x,datatype);
        else if (i == 1) image_sptr->write(filename_y,datatype);
        else if (i == 2) image_sptr->write(filename_z,datatype);
    }
}

template<class dataType>
void NiftiImageData3DTensor<dataType>::flip_component(const int dim)
{
    std::cout << "\nFlipping multicomponent image in dim number: " << dim << "..." << std::flush;

    // Check the dimension to flip, that dims==5 and nu==3
    if (dim < 0 || dim > 2)
        throw std::runtime_error("\n\tDimension to flip should be between 0 and 2.");

    // Data is ordered such that the multicomponent info is last.
    // So, the first third of the data is the x-values, second third is y and last third is z.
    // Start index is therefore = dim_number * num_voxels/3
    int start_index =   dim   * int(this->_nifti_image->nvox/3);
    int end_index   = (dim+1) * int(this->_nifti_image->nvox/3 - 1);

    for (int i=start_index; i<=end_index; i++)
        this->_data[i] = -this->_data[i];
}

namespace sirf {
template class NiftiImageData3DTensor<float>;
}
