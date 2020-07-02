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
\brief Class for deformation field transformations.

\author Richard Brown
\author SyneRBI
*/

#include "sirf/Reg/NiftiImageData3DBSpline.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/NiftyMoMo/BSplineTransformation.h"

using namespace sirf;

template<class dataType>
NiftiImageData3DBSpline<dataType>::NiftiImageData3DBSpline(const NiftiImageData3DDeformation<dataType> &def, float spacing[])
{
    // Get any of the tensor components as a 3d image
    nifti_image *ref_ptr = def.get_tensor_component(0)->get_raw_nifti_sptr().get();
    // Create the NiftyMoMo bspline transformation class
    NiftyMoMo::BSplineTransformation bspline(ref_ptr, 1, spacing);
    // Convert DVF to CPG
    bspline.GetDVFGradientWRTTransformationParameters(def.clone()->get_raw_nifti_sptr().get(), ref_ptr);
    // Get output
    nifti_image *cpg_ptr = bspline.GetTransformationAsImage();
    cpg_ptr->intent_p1 = SPLINE_VEL_GRID;
    *this = NiftiImageData3DBSpline<dataType>(*cpg_ptr);
    this->check_dimensions(NiftiImageData<dataType>::_3DBSpl);
}

template<class dataType>
void NiftiImageData3DBSpline<dataType>::create_from_3D_image(const NiftiImageData<dataType> &image)
{
    NiftiImageData3DTensor<dataType>::create_from_3D_image(image);
    this->_nifti_image->intent_p1 = SPLINE_VEL_GRID;
}

template<class dataType>
NiftiImageData3DDeformation<dataType> NiftiImageData3DBSpline<dataType>::get_as_deformation_field(const NiftiImageData<dataType> &ref) const
{
    // Get spacing of reference image
    float spacing[3];
    for (unsigned i=0; i<3; ++i)
        spacing[i] = this->_nifti_image->pixdim[i+1];
    // Create the NiftyMoMo bspline transformation class
    NiftyMoMo::BSplineTransformation bspline(ref.clone()->get_raw_nifti_sptr().get(), 1, spacing);
    // Set the CPG
    bspline.SetParameters(static_cast<float*>(this->_nifti_image->data), false);
    // Get the DVF
    nifti_image *output_def_ptr = bspline.GetDeformationVectorField(ref.get_raw_nifti_sptr().get());
    return NiftiImageData3DDeformation<dataType>(*output_def_ptr);
}

template<class dataType>
NiftiImageData3DBSpline<dataType>*
NiftiImageData3DBSpline<dataType>::get_inverse_impl_nr(const std::shared_ptr<const NiftiImageData<dataType> >) const
{
    throw std::runtime_error("NiftiImageData3DBSpline::get_inverse_impl_nr not yet implemented.");
}

template<class dataType>
NiftiImageData3DBSpline<dataType>*
NiftiImageData3DBSpline<dataType>::get_inverse_impl_vtk(const std::shared_ptr<const NiftiImageData<dataType> >) const
{
    throw std::runtime_error("NiftiImageData3DBSpline::get_inverse_impl_vtk not yet implemented.");
#ifndef SIRF_VTK
    throw std::runtime_error("Build SIRF with VTK support for this functionality");
#endif
}

namespace sirf {
template class NiftiImageData3DBSpline<float>;
}
