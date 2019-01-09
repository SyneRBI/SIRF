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

#ifndef _NIFTIIMAGEDATA3DDEFORMATION_H_
#define _NIFTIIMAGEDATA3DDEFORMATION_H_

#include "sirf/cReg/NiftiImageData3DTensor.h"
#include "sirf/cReg/SIRFRegTransformation.h"

namespace sirf {

// Forward declarations
template<class dataType> class NiftiImageData3D;
template<class dataType> class NiftiImageData3DDisplacement;

/// SIRF nifti image data deformation field image
template<class dataType>
class NiftiImageData3DDeformation : public NiftiImageData3DTensor<dataType>, public SIRFRegTransformation<dataType>
{
public:
    /// Constructor
    NiftiImageData3DDeformation() {}

    /// Filename constructor
    NiftiImageData3DDeformation(const std::string &filename)
        : NiftiImageData3DTensor<dataType>(filename) { this->check_dimensions(this->_3DDef); }

    /// Nifti constructor
    NiftiImageData3DDeformation(const nifti_image &image_nifti)
        : NiftiImageData3DTensor<dataType>(image_nifti) { this->check_dimensions(this->_3DDef); }

    /// Construct from general tensor
    NiftiImageData3DDeformation(const NiftiImageData<dataType>& tensor)
        : NiftiImageData3DTensor<dataType>(tensor) { this->check_dimensions(this->_3DDef); }

    /// Create from 3 individual components
    NiftiImageData3DDeformation(const NiftiImageData3D<dataType> &x, const NiftiImageData3D<dataType> &y, const NiftiImageData3D<dataType> &z)
        : NiftiImageData3DTensor<dataType>(x,y,z) { this->check_dimensions(this->_3DDef); }

    /// Create from displacement field image
    NiftiImageData3DDeformation(const NiftiImageData3DDisplacement<dataType> &disp);

    /// Create from 3D image
    void create_from_3D_image(const NiftiImageData3D<dataType> &image);

    /// Create from CPP image
    void create_from_cpp(NiftiImageData3DTensor<dataType> &cpp, const NiftiImageData3D<dataType> &ref);

    /// Get as deformation field
    virtual NiftiImageData3DDeformation get_as_deformation_field(const NiftiImageData3D<dataType> &ref) const;

    /// Compose multiple transformations into single deformation field
    static NiftiImageData3DDeformation compose_single_deformation(const std::vector<const SIRFRegTransformation<dataType> *> &transformations, const NiftiImageData3D<dataType> &ref);

    /// Compose multiple transformations into single deformation field
    static NiftiImageData3DDeformation compose_single_deformation(const std::vector<std::shared_ptr<const SIRFRegTransformation<dataType> > > &transformations, const NiftiImageData3D<dataType> &ref);

    virtual NiftiImageData3DDeformation* same_image_data() const
    {
        return new NiftiImageData3DDeformation;
    }
    /// Write
    virtual void write(const std::string &filename) const { this->NiftiImageData<dataType>::write(filename); }
    /// Clone and return as unique pointer.
    std::unique_ptr<NiftiImageData3DDeformation> clone() const
    {
	return std::unique_ptr<NiftiImageData3DDeformation>(this->clone_impl());
    }
protected:
    /// Clone helper function. Don't use.
    virtual NiftiImageData3DDeformation* clone_impl() const
    {
	return new NiftiImageData3DDeformation(*this);
    }
};
}

#endif
