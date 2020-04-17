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
\brief Resampling class based on nifty resample

\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include <nifti1_io.h>
#include <string>
#include <vector>
#include <iostream>
#include "sirf/Reg/Resample.h"
#include "sirf/iUtilities/iutilities.h"

namespace NiftyMoMo {
class BSplineTransformation;
}

namespace sirf {

namespace detail {
/*! \ingroup Registration
  \brief This is an internal class requied by NiftyResample to handle complex images.

NiftyReg doesn't currently handle complex nifti images. So we split complex (e.g., MR)
images into two images, a real and imaginary component. We can then resample them and
then join the two parts back together.
 */
template<class dataType>
class ComplexNiftiImageData
{
public:
    /// Constructor
    ComplexNiftiImageData() {}
    /// Destructor
    virtual ~ComplexNiftiImageData() {}
    /// is complex
    bool is_complex() const { return _imag_sptr != nullptr; }
    /// Get real component
    const std::shared_ptr<const NiftiImageData<dataType> > real() const { return _real_sptr; }
    /// Get real component
    std::shared_ptr<NiftiImageData<dataType> > &real() { return _real_sptr; }
    /// Get imaginary component
    const std::shared_ptr<const NiftiImageData<dataType> > imag() const { return _imag_sptr; }
    /// Get imaginary component
    std::shared_ptr<NiftiImageData<dataType> > &imag() { return _imag_sptr; }
    /// size
    size_t size() const { return is_complex() ? 2 : 1; }
    /// at
    const std::shared_ptr<const NiftiImageData<dataType> > at(const unsigned idx) const { check_bounds(idx); return idx == 0 ? real() : imag(); }
    /// at
    std::shared_ptr<NiftiImageData<dataType> > at(const unsigned idx) { check_bounds(idx); return idx == 0 ? real() : imag(); }

private:
    void check_bounds(const unsigned idx) const
    {
        if (idx > 1)
            throw std::runtime_error("ComplexNiftiImageData::at(): Exceeds index range");
        if (idx == 1 && !is_complex())
            throw std::runtime_error("ComplexNiftiImageData::at(): Trying to access imaginary part of non-complex image");
    }
    std::shared_ptr<NiftiImageData<dataType> > _real_sptr;
    std::shared_ptr<NiftiImageData<dataType> > _imag_sptr;
};
}

/*!
\ingroup Registration
\brief Resampling class based on nifty resample

The reference image and floating image can have nt and/or nu != 1.

\author Richard Brown
\author CCP PETMR
*/
template<class dataType>
class NiftyResample : public Resample<dataType>
{
public:

    /// Constructor
    NiftyResample() {}

    /// Destructor
    virtual ~NiftyResample() {}

    /// Process
    DEPRECATED virtual void process();

    /// Do the forward transformation
    virtual std::shared_ptr<ImageData> forward(const std::shared_ptr<const ImageData> input_sptr);

    /// Do the forward transformation
    virtual void forward(std::shared_ptr<ImageData> output_sptr, const std::shared_ptr<const ImageData> input_sptr);

    /// Do the adjoint transformation
    virtual std::shared_ptr<ImageData> adjoint(const std::shared_ptr<const ImageData> input_sptr);

    /// Do the adjoint transformation
    virtual void adjoint(std::shared_ptr<ImageData> output_sptr, const std::shared_ptr<const ImageData> input_sptr);

protected:

    /// Set up
    virtual void set_up();

    /// Set up forward
    virtual void set_up_forward();

    /// Set up adjoint
    virtual void set_up_adjoint();

    /// Set up the input images (convert from ImageData to NiftiImageData if necessary)
    void set_up_input_images();

    /// Reference image as a NiftiImageData
    detail::ComplexNiftiImageData<dataType> _reference_image_niftis;
    /// Floating image as a NiftiImageData
    detail::ComplexNiftiImageData<dataType> _floating_image_niftis;
    /// Forward resampled image as a NiftiImageData
    detail::ComplexNiftiImageData<dataType> _output_image_forward_niftis;
    /// Adjoint resampled image as a NiftiImageData
    detail::ComplexNiftiImageData<dataType> _output_image_adjoint_niftis;

    /// Deformation
    std::shared_ptr<NiftiImageData3DDeformation<dataType> > _deformation_sptr;
    /// Needed for the adjoint transformation
    std::shared_ptr<NiftyMoMo::BSplineTransformation> _adjoint_transformer_sptr;
    /// Adjoint reference weights. Vector as may be complex
    std::shared_ptr<NiftiImageData<dataType> > _adjoint_input_weights_sptr;
    /// Adjoint output weights. Vector as may be complex
    std::shared_ptr<NiftiImageData<dataType> > _adjoint_output_weights_sptr;
};
}
