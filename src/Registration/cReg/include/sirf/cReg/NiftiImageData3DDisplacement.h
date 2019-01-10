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
\brief Class for deformation/displacement SIRF image data.

\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include "sirf/cReg/NiftiImageData3DTensor.h"
#include "sirf/cReg/NiftiImageData3DDeformation.h"
#include "sirf/cReg/Transformation.h"

namespace sirf {
template<class dataType> class NiftiImageData3D;

/// SIRF nifti image data displacement field image
template<class dataType>
class NiftiImageData3DDisplacement : public NiftiImageData3DTensor<dataType>, public Transformation<dataType>
{
public:
    /// Constructor
    NiftiImageData3DDisplacement() {}

    /// Filename constructor
    NiftiImageData3DDisplacement(const std::string &filename)
        : NiftiImageData3DTensor<dataType>(filename) { this->check_dimensions(this->_3DDisp); }

    /// Nifti constructor
    NiftiImageData3DDisplacement(const nifti_image &image_nifti)
        : NiftiImageData3DTensor<dataType>(image_nifti) { this->check_dimensions(this->_3DDisp); }

    /// Construct from general tensor
    NiftiImageData3DDisplacement(const NiftiImageData<dataType>& tensor)
        : NiftiImageData3DTensor<dataType>(tensor) { this->check_dimensions(this->_3DDisp); }

    /// Create from 3 individual components
    NiftiImageData3DDisplacement(const NiftiImageData3D<dataType> &x, const NiftiImageData3D<dataType> &y, const NiftiImageData3D<dataType> &z)
        : NiftiImageData3DTensor<dataType>(x,y,z) { this->_nifti_image->intent_p1 = 1; }

    /// Create from deformation field image
    NiftiImageData3DDisplacement(const NiftiImageData3DDeformation<dataType> &def);

    /// Create from 3D image
    void create_from_3D_image(const NiftiImageData3D<dataType> &image);

    /// Get as deformation field
    virtual NiftiImageData3DDeformation<dataType> get_as_deformation_field(const NiftiImageData3D<dataType> &ref) const;

    virtual NiftiImageData3DDisplacement* same_image_data() const
    {
        return new NiftiImageData3DDisplacement;
    }
    /// Write
    virtual void write(const std::string &filename) const { this->NiftiImageData<dataType>::write(filename); }
    /// Clone and return as unique pointer.
    std::unique_ptr<NiftiImageData3DDisplacement> clone() const
    {
	return std::unique_ptr<NiftiImageData3DDisplacement>(this->clone_impl());
    }
protected:
    /// Clone helper function. Don't use.
    virtual NiftiImageData3DDisplacement* clone_impl() const
    {
	return new NiftiImageData3DDisplacement<dataType>(*this);
    }
};
}
