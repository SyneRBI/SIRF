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
\brief Class for getting image gradient WRT a transformation and multiplying by image.

\author Richard Brown
\author SyneRBI
*/

#pragma once

#include <memory>

namespace sirf {

template<class dataType> class NiftyResample;
template<class dataType> class NiftiImageData3DDeformation;
class ImageData;

/*!
\ingroup Registration
\brief Class for getting image gradient WRT a transformation and multiplying by image.

\author Richard Brown
\author SyneRBI
*/
template<class dataType>
class ImageGradientWRTDeformationTimesImage
{
public:

    /// Constructor
    ImageGradientWRTDeformationTimesImage();

    /// Set the resampler
    void set_resampler(const std::shared_ptr<NiftyResample<dataType> > resampler_sptr);

    /// Forward in place (resample image)
    virtual void forward(std::shared_ptr<ImageData> &im_sptr, const std::shared_ptr<const NiftiImageData3DDeformation<dataType> > &deformation_sptr);

    /// Forward (resample image)
    virtual std::shared_ptr<const ImageData> forward(const std::shared_ptr<const NiftiImageData3DDeformation<dataType> > deformation_sptr);

    /// Backward in place (get image gradient wrt transformation)
    virtual void backward(std::shared_ptr<NiftiImageData3DDeformation<dataType> > &output_transformation_sptr, const std::shared_ptr<const ImageData> image_to_multiply_sptr);

    /// Backward (get image gradient wrt transformation)
    virtual std::shared_ptr<const NiftiImageData3DDeformation<dataType> > backward(const std::shared_ptr<const ImageData> image_to_multiply_sptr);

private:

    /// Resampler
    std::shared_ptr<NiftyResample<dataType> > _resampler_sptr;
};
}
