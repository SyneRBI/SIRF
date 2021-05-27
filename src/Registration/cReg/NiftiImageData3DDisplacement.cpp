/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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
\author SyneRBI
*/

#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include <_reg_localTrans.h>
#ifdef SIRF_VTK
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkImageBSplineCoefficients.h"
#include "vtkBSplineTransform.h"
#include "vtkTransformToGrid.h"
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

template<class dataType>
NiftiImageData3DDisplacement<dataType>*
NiftiImageData3DDisplacement<dataType>::get_inverse_impl_nr(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr) const
{
    // Convert to deformation
    NiftiImageData3DDeformation<dataType> def(*this);

    // Inverse with NR
    auto def_inv = def.get_inverse(image_sptr, false);

    // Convert back to disp
    return new NiftiImageData3DDisplacement<dataType>(*def_inv);
}

template<class dataType>
NiftiImageData3DDisplacement<dataType>*
NiftiImageData3DDisplacement<dataType>::get_inverse_impl_vtk(const std::shared_ptr<const NiftiImageData<dataType> > image_sptr) const
{
#ifndef SIRF_VTK
    throw std::runtime_error("Build SIRF with VTK support for this functionality");
#else
    if (image_sptr != nullptr)
        throw std::runtime_error("VTK inversion not yet implemented for change of grid. You should simply be able to resample.");
    if (typeid(dataType) != typeid(float))
        throw std::runtime_error("VTK inversion not yet implemented for dataTypes other than float.");

    // Get image info
    const int *dims = this->get_dimensions();
    const float *spacing = this->get_raw_nifti_sptr()->pixdim;

    // Initialise a vtkImageData
    auto image_data_sptr = vtkSmartPointer<vtkImageData>::New();
    image_data_sptr->SetDimensions(dims[1],dims[2],dims[3]);
    image_data_sptr->SetSpacing(spacing[1],spacing[2],spacing[3]);
    image_data_sptr->AllocateScalars(VTK_FLOAT,dims[5]);

    // Copy NiftiImageData to vtkImageData
    int idx[7] = {0, 0, 0, 0, 0, 0, 0};
    for (idx[0]=0; idx[0]<dims[1]; ++idx[0])
        for (idx[1]=0; idx[1]<dims[2]; ++idx[1])
            for (idx[2]=0; idx[2]<dims[3]; ++idx[2])
                for (idx[4]=0; idx[4]<dims[5]; ++idx[4]) // t is 3, skip to 4 for tensor component
                    image_data_sptr->SetScalarComponentFromFloat(
                                idx[0],idx[1],idx[2],idx[4],(*this)(idx));

    // create b-spline coefficients for the displacement grid image_data_sptr
    auto bspline_filter = vtkSmartPointer<vtkImageBSplineCoefficients>::New();
    bspline_filter->SetInputData(image_data_sptr);
    bspline_filter->Update();

    // use these b-spline coefficients to create a transform
    auto bspline_transform = vtkSmartPointer<vtkBSplineTransform>::New();
    bspline_transform->SetCoefficientData(bspline_filter->GetOutput());

    // invert the b-spline transform onto a new grid
    auto grid_maker = vtkSmartPointer<vtkTransformToGrid>::New();
    grid_maker->SetInput(bspline_transform->GetInverse());
    grid_maker->SetGridOrigin(image_data_sptr->GetOrigin());
    grid_maker->SetGridSpacing(image_data_sptr->GetSpacing());
    grid_maker->SetGridExtent(image_data_sptr->GetExtent());
    grid_maker->SetGridScalarTypeToFloat();
    grid_maker->Update();

    // Get inverse displacement as an image
    auto inverse_image = grid_maker->GetOutput();

    // Copy vtkImageData back to NiftiImageData
    auto *output_ptr = new NiftiImageData3DDisplacement<dataType>(*this);
    for (idx[0]=0; idx[0]<dims[1]; ++idx[0])
        for (idx[1]=0; idx[1]<dims[2]; ++idx[1])
            for (idx[2]=0; idx[2]<dims[3]; ++idx[2])
                for (idx[4]=0; idx[4]<dims[5]; ++idx[4])// t is 3, skip to 4 for tensor component
                    (*output_ptr)(idx) = inverse_image->GetScalarComponentAsFloat(idx[0],idx[1],idx[2],idx[4]);

    return output_ptr;
#endif
}

namespace sirf {
template class NiftiImageData3DDisplacement<float>;
}
