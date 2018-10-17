/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2018 University College London

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

#ifndef SIRFREG_PARAMETERS_HANDLERS
#define SIRFREG_PARAMETERS_HANDLERS

#include "data_handle.h"

namespace sirf {

    // ---------------------------------------------------------------------------- //
    // NiftiImageData
    // ---------------------------------------------------------------------------- //
	void*
        cSIRFReg_NiftiImageDataParameter(const DataHandle* handle, const char* name);

    // ---------------------------------------------------------------------------- //
    // SIRFReg
    // ---------------------------------------------------------------------------- //
	void*
		cSIRFReg_setSIRFRegParameter(void* hp, const char* name, const void* hv);
	void*
		cSIRFReg_SIRFRegParameter(const DataHandle* handle, const char* name);

    // ---------------------------------------------------------------------------- //
    // SIRFRegNiftyF3dSym
    // ---------------------------------------------------------------------------- //
    void*
        cSIRFReg_setSIRFRegNiftyF3dSymParameter(void* hp, const char* name, const void* hv);

    // ---------------------------------------------------------------------------- //
    // SIRFRegNiftyResample
    // ---------------------------------------------------------------------------- //
    void*
        cSIRFReg_setSIRFRegNiftyResampleParameter(void* hp, const char* name, const void* hv);
    void*
        cSIRFReg_SIRFRegNiftyResampleParameter(const DataHandle* handle, const char* name);

    // ---------------------------------------------------------------------------- //
    // SIRFRegImageWeightedMean
    // ---------------------------------------------------------------------------- //
    void*
        cSIRFReg_SIRFRegImageWeightedMeanParameter(const DataHandle* handle, const char* name);

    // ---------------------------------------------------------------------------- //
    // SIRFRegMat44
    // ---------------------------------------------------------------------------- //
    void*
        cSIRFReg_SIRFRegMat44Parameter(const DataHandle* handle, const char* name);
}

#endif
