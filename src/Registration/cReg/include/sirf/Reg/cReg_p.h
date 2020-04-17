/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2017 - 2020 University College London

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

#pragma once

#include "sirf/iUtilities/DataHandle.h"

namespace sirf {

    // ---------------------------------------------------------------------------- //
    // NiftiImageData
    // ---------------------------------------------------------------------------- //
    void* cReg_NiftiImageDataParameter(const DataHandle* handle, const char* name);

    // ---------------------------------------------------------------------------- //
    // Registration
    // ---------------------------------------------------------------------------- //
    void* cReg_setRegistrationParameter(void* hp, const char* name, const void* hv);

    // ---------------------------------------------------------------------------- //
    // NiftyRegistration
    // ---------------------------------------------------------------------------- //
    void* cReg_setNiftyRegistrationParameter(void* hp, const char* name, const void* hv);

    // ---------------------------------------------------------------------------- //
    // NiftyF3dSym
    // ---------------------------------------------------------------------------- //
    void* cReg_setNiftyF3dSymParameter(void* hp, const char* name, const void* hv);

#ifdef SIRF_SPM
    // ---------------------------------------------------------------------------- //
    // SPMRegistration
    // ---------------------------------------------------------------------------- //
    void* cReg_setSPMRegistrationParameter(void* hp, const char* name, const void* hv);
#endif

    // ---------------------------------------------------------------------------- //
    // NiftyResample
    // ---------------------------------------------------------------------------- //
    void* cReg_setNiftyResampleParameter(void* hp, const char* name, const void* hv);
    void* cReg_NiftyResampleParameter(const DataHandle* handle, const char* name);

    // ---------------------------------------------------------------------------- //
    // ImageWeightedMean
    // ---------------------------------------------------------------------------- //
    void* cReg_ImageWeightedMeanParameter(const DataHandle* handle, const char* name);

    // ---------------------------------------------------------------------------- //
    // AffineTransformation
    // ---------------------------------------------------------------------------- //
    void* cReg_AffineTransformationParameter(const DataHandle* handle, const char* name);
}
