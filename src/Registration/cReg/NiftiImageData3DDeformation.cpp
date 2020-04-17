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
\brief Class for deformation field transformations.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/Transformation.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include <_reg_globalTrans.h>
#include <sstream>
#include <_reg_localTrans.h>

using namespace sirf;

template<class dataType>
NiftiImageData3DDeformation<dataType>::NiftiImageData3DDeformation(const NiftiImageData3DDisplacement<dataType> &disp)
{
    // Get the def field from the disp field
    NiftiImageData3DTensor<dataType> temp = disp;
    reg_getDeformationFromDisplacement(temp.get_raw_nifti_sptr().get());
    temp.get_raw_nifti_sptr()->intent_p1 = DEF_FIELD;
    *this = temp;
}

template<class dataType>
void NiftiImageData3DDeformation<dataType>::create_from_3D_image(const NiftiImageData<dataType> &image)
{
    this->NiftiImageData3DTensor<dataType>::create_from_3D_image(image);
    //_nifti_image->intent_p1 = 0; not necessary. 0 by default
}

template<class dataType>
void NiftiImageData3DDeformation<dataType>::create_from_cpp(NiftiImageData3DTensor<dataType> &cpp, const NiftiImageData<dataType> &ref)
{
    this->create_from_3D_image(ref);

    reg_spline_getDeformationField(cpp.get_raw_nifti_sptr().get(),
                                   this->_nifti_image.get(),
                                   NULL,
                                   false, //composition
                                   true // bspline
                                   );
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DDeformation<dataType>::get_as_deformation_field(const NiftiImageData<dataType> &ref) const
{
    NiftiImageData3DDeformation<dataType> output_def;
    output_def.create_from_3D_image(ref);
    nifti_image * def_ptr = output_def.get_raw_nifti_sptr().get();

    // Initialise the deformation field with an identity transformation
    reg_tools_multiplyValueToImage(def_ptr,def_ptr,0.f);
    reg_getDeformationFromDisplacement(def_ptr);
    def_ptr->intent_p1=DEF_FIELD;

    // Not marked const so have to copy unfortunately
    std::shared_ptr<NiftiImageData3DDeformation<dataType> > copy_of_input_def_sptr =
            this->clone();

    reg_defField_compose(copy_of_input_def_sptr->get_raw_nifti_sptr().get(),
                         def_ptr,
                         nullptr);
    return output_def;
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DDeformation<dataType>::compose_single_deformation(const std::vector<const Transformation<dataType>*> &transformations, const NiftiImageData<dataType> &ref)
{
    if (transformations.size() == 0)
        throw std::runtime_error("NiftiImageData3DDeformation::compose_single_deformation no transformations given.");

    NiftiImageData3DDeformation def = transformations.at(0)->get_as_deformation_field(ref);

    for (unsigned i=1; i<transformations.size(); ++i) {
        NiftiImageData3DDeformation temp = transformations.at(i)->get_as_deformation_field(ref);
        reg_defField_compose(temp.get_raw_nifti_sptr().get(),def.get_raw_nifti_sptr().get(),nullptr);
    }
    return def;
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DDeformation<dataType>::compose_single_deformation(const std::vector<std::shared_ptr<const Transformation<dataType> > > &transformations, const NiftiImageData<dataType> &ref)
{
    std::vector<const Transformation<dataType>*> vec;
    for (unsigned i=0; i<transformations.size(); ++i)
        vec.push_back(transformations.at(i).get());
    return compose_single_deformation(vec, ref);
}

template<class dataType>
std::shared_ptr<const NiftiImageData3DDeformation<dataType> > NiftiImageData3DDeformation<dataType>::get_inverse(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr) const
{
    // Output is based on original deformation field or some other
    // input, depending on whether input argument is nullptr or not
    nifti_image *output_ptr;
    if (image_sptr == nullptr)
        output_ptr=nifti_copy_nim_info(this->get_raw_nifti_sptr().get());
    else
        output_ptr=nifti_copy_nim_info(image_sptr->get_raw_nifti_sptr().get());
    output_ptr->ndim = output_ptr->dim[0] = 5;
    output_ptr->nt = output_ptr->dim[4] = 1;
    output_ptr->nu = output_ptr->dim[5] = output_ptr->nz>1 ? 3 : 2;
    output_ptr->nvox = size_t(output_ptr->nx *
       output_ptr->ny * output_ptr->nz *
       output_ptr->nt * output_ptr->nu);
    output_ptr->nbyper = this->get_raw_nifti_sptr()->nbyper;
    output_ptr->datatype = this->get_raw_nifti_sptr()->datatype;
    output_ptr->intent_code = NIFTI_INTENT_VECTOR;
    output_ptr->intent_p1 = this->get_raw_nifti_sptr()->intent_p1;
    output_ptr->intent_p2 = this->get_raw_nifti_sptr()->intent_p2;
    output_ptr->scl_slope = 1.f;
    output_ptr->scl_inter = 0.f;
    output_ptr->data = static_cast<void *>(malloc
       (output_ptr->nvox*size_t(output_ptr->nbyper)));
    // Method not marked const, so have to clone the deformation field
    reg_defFieldInvert(this->clone()->get_raw_nifti_sptr().get(),output_ptr,1.0e-6f);

    // Create shared_ptr, delete the original ptr
    std::shared_ptr<const NiftiImageData3DDeformation<dataType> > output_sptr =
            std::make_shared<NiftiImageData3DDeformation<dataType> >(*output_ptr);
    nifti_image_free(output_ptr);
    return output_sptr;
}

namespace sirf {
template class NiftiImageData3DDeformation<float>;
}
