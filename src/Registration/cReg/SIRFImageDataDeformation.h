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

#ifndef _SIRFIMAGEDATADEFORMATION_H_
#define _SIRFIMAGEDATADEFORMATION_H_


#include "SIRFImageData.h"

/// SIRF image data
class SIRFImageDataDeformation : public SIRFImageData
{
public:

    /// Constructor
    SIRFImageDataDeformation() {}

    /// Filename constructor
    SIRFImageDataDeformation(const std::string &filename);

    /// Nifti constructor
    SIRFImageDataDeformation(const nifti_image *image_nifti);

    /// Nifti shared_ptr constructor
    SIRFImageDataDeformation(const std::shared_ptr<nifti_image> image_nifti);

    /// Create from 3D image
    void create_from_3D_image(const SIRFImageData &image);

    /// Save to file
    void save_to_file(const std::string &filename, bool split_xyz=false, std::string type = "") const;

    /// Deep copy
    SIRFImageDataDeformation deep_copy() const;

protected:


};

#endif
