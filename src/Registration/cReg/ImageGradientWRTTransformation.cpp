/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 University College London

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
\brief Class for getting image gradient WRT a transformation.

\author Richard Brown
\author SyneRBI
*/

#include "sirf/Reg/ImageGradientWRTTransformation.h"
#include "sirf/Reg/Resample.h"

using namespace sirf;

template<class dataType>
void
ImageGradientWRTTransformation<dataType>::
set_resampler(const std::shared_ptr<Resample<dataType> > resampler_sptr)
{
    _resampler_sptr = resampler_sptr;
}

template<class dataType>
void
ImageGradientWRTTransformation<dataType>::
forward(std::shared_ptr<Transformation<dataType> > &output_transformation_sptr, const std::shared_ptr<const ImageData> source_im_sptr)
{
    _resampler_sptr->get_image_gradient_wrt_transformation(output_transformation_sptr, source_im_sptr);
}

template<class dataType>
std::shared_ptr<const Transformation<dataType> >
ImageGradientWRTTransformation<dataType>::
forward(const std::shared_ptr<const ImageData> source_im_sptr)
{
    return std::move(_resampler_sptr->get_image_gradient_wrt_transformation(source_im_sptr));
}

namespace sirf {
template class ImageGradientWRTTransformation<float>;
}
