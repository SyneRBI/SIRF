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

#pragma once

#include <memory>

namespace sirf {

template<class dataType> class Resample;
template<class dataType> class Transformation;
class ImageData;

/*!
\ingroup Registration
\brief Class for converting control point grids to deformation field transformations.

\author Richard Brown
\author SyneRBI
*/
template<class dataType>
class ImageGradientWRTTransformation
{
public:

    /// Constructor
    ImageGradientWRTTransformation();

    /// Set the resampler
    void set_resampler(const std::shared_ptr<Resample<dataType> > resampler_sptr);

    /// Forward in place (get image gradient wrt transformation)
    virtual void forward(std::shared_ptr<Transformation<dataType> > &output_transformation_sptr, const std::shared_ptr<const ImageData> source_im_sptr);

    /// Forward (get image gradient wrt transformation)
    virtual std::shared_ptr<const Transformation<dataType> > forward(const std::shared_ptr<const ImageData> source_im_sptr);

private:

    /// Resampler
    std::shared_ptr<Resample<dataType> > _resampler_sptr;
};
}
