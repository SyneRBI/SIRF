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
\brief Base class for SIRF image data.

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFIMAGEDATA_H_
#define _SIRFIMAGEDATA_H_

#include <memory>
#include <string>
#include <nifti1_io.h>
#include <stir_data_containers.h>

/// SIRF image data
class SIRFImageData
{
public:

    /// Constructor
    SIRFImageData() {}

    /// Destructor
    virtual ~SIRFImageData() {}

    /// Set image as nifti
    void set_image_nifti(const nifti_image &image_nifti);

    /// Set image as filename
    void set_image_filename(const std::string &image_filename);

    /// Set image as SIRF's PETImageData
    void set_image_PETImageData(const PETImageData &image_sirf_pet);

    /// Get image as nifti
    nifti_image *get_image_as_nifti();

protected:

    /// Image data as a nifti object
    std::shared_ptr<nifti_image>  _image_nifti;
    /// Image data as a filename
    std::string                   _image_filename;
    /// Image data as a SIRF Image Data
    std::shared_ptr<PETImageData> _image_sirf_pet;
};

#endif
