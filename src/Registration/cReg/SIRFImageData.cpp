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

#include "SIRFImageData.h"
#include "SIRFRegMisc.h"

using namespace std;

void SIRFImageData::set_image_nifti(const nifti_image &image_nifti)
{
    _image_nifti = make_shared<nifti_image>(image_nifti);
    _image_filename = "";
    _image_sirf_pet.reset();
}

void SIRFImageData::set_image_filename(const std::string &image_filename)
{
    _image_filename = image_filename;
    _image_nifti.reset();
    _image_sirf_pet.reset();
}

void SIRFImageData::set_image_PETImageData(const PETImageData &image_sirf_pet)
{
    _image_sirf_pet = make_shared<PETImageData>(image_sirf_pet);
    _image_nifti.reset();
    _image_filename = "";
}

nifti_image *SIRFImageData::get_image_as_nifti()
{
    // If the image is nifti
    if (_image_nifti) return _image_nifti.get();

    // If the image is filename
    else /*if (_image_filename != "") */{
        SIRFRegMisc::open_nifti_image(_image_nifti, _image_filename);
        return _image_nifti.get();
    }

}
