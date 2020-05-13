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
\brief Class for displacement field transformations.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include <_reg_localTrans.h>
#ifdef SIRF_VTK
#include "vtkSmartPointer.h"
#include "vtkGridTransform.h"
#include "vtkImageData.h"
#endif

using namespace sirf;

template<class dataType>
NiftiImageData3DDisplacement<dataType>::NiftiImageData3DDisplacement(const NiftiImageData3DDeformation<dataType> &def)
{
    // Get the disp field from the def field
    NiftiImageData3DTensor<dataType> temp = def;
    reg_getDisplacementFromDeformation(temp.get_raw_nifti_sptr().get());
    temp.get_raw_nifti_sptr()->intent_p1 = DISP_FIELD;
    *this = temp;
}

template<class dataType>
void NiftiImageData3DDisplacement<dataType>::create_from_3D_image(const NiftiImageData<dataType> &image)
{
    this->NiftiImageData3DTensor<dataType>::create_from_3D_image(image);
    this->_nifti_image->intent_p1 = 1;
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DDisplacement<dataType>::get_as_deformation_field(const NiftiImageData<dataType> &ref) const
{
    NiftiImageData3DDeformation<dataType> output_def;
    output_def.create_from_3D_image(ref);
    nifti_image * def_ptr = output_def.get_raw_nifti_sptr().get();

    // Initialise the deformation field with an identity transformation
    reg_tools_multiplyValueToImage(def_ptr,def_ptr,0.f);
    reg_getDeformationFromDisplacement(def_ptr);
    def_ptr->intent_p1=DEF_FIELD;

    // Not marked const so have to copy unfortunately
    std::shared_ptr<NiftiImageData3DDisplacement<dataType> > copy_of_input_disp_sptr =
            this->clone();
    // Convert displacement to deformation
    reg_getDeformationFromDisplacement(copy_of_input_disp_sptr->get_raw_nifti_sptr().get());
    reg_defField_compose(copy_of_input_disp_sptr->get_raw_nifti_sptr().get(),
                         def_ptr,
                         nullptr);
    return output_def;
}

#ifdef SIRF_VTK
template<class dataType>
NiftiImageData3DDisplacement<dataType>
NiftiImageData3DDisplacement<dataType>::get_inverse_vtk() const
{
    const int *dims = this->get_dimensions();
    const float *spacing = this->get_raw_nifti_sptr()->pixdim;

    using vtkImage = vtkSmartPointer<vtkImageData>;
    using vtkTransform = vtkSmartPointer<vtkGridTransform>;

    vtkImage image_data_sptr = vtkImage::New();
    image_data_sptr->SetDimensions(dims[1],dims[2],dims[3]);
    image_data_sptr->SetSpacing(spacing[1],spacing[2],spacing[3]);
    image_data_sptr->AllocateScalars(VTK_FLOAT,dims[5]);

    int idx[7] = {0, 0, 0, 0, 0, 0, 0};
    for (idx[0]=0; idx[0]<dims[1]; ++idx[0])
        for (idx[1]=0; idx[1]<dims[2]; ++idx[1])
            for (idx[2]=0; idx[2]<dims[3]; ++idx[2])
                for (idx[4]=0; idx[4]<dims[5]; ++idx[4]) // t is 3, skip to 4 for tensor component
                    image_data_sptr->SetScalarComponentFromFloat(
                                idx[0],idx[1],idx[2],idx[4],(*this)(idx));

    vtkTransform transform_sptr = vtkTransform::New();
    transform_sptr->SetDisplacementGridData(image_data_sptr);
    transform_sptr->Update();
    vtkTransform inverse_transform_sptr =
            vtkGridTransform::SafeDownCast(transform_sptr->GetInverse());

    vtkImage inverse_image = inverse_transform_sptr->GetDisplacementGrid();

    NiftiImageData3DDisplacement<float> output = *this->clone();
    for (idx[0]=0; idx[0]<dims[1]; ++idx[0])
        for (idx[1]=0; idx[1]<dims[2]; ++idx[1])
            for (idx[2]=0; idx[2]<dims[3]; ++idx[2])
                for (idx[4]=0; idx[4]<dims[5]; ++idx[4]) {// t is 3, skip to 4 for tensor component
                    output(idx) = inverse_image->GetScalarComponentAsFloat(idx[0],idx[1],idx[2],idx[4]);
                }

    return output;
}
#endif

namespace sirf {
template class NiftiImageData3DDisplacement<float>;
}
