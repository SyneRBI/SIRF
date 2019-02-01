/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2017 University College London.
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
#ifndef CSIRF_TO_MATLAB_INTERFACE
#define CSIRF_TO_MATLAB_INTERFACE

#define CSIRF_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

#ifndef CSIRF_FOR_MATLAB
#define PTR_INT size_t
#define PTR_FLOAT size_t
#define PTR_DOUBLE size_t
 extern "C" {
#else
#define PTR_INT int*
#define PTR_FLOAT float*
#define PTR_DOUBLE double*
#endif
EXPORTED_FUNCTION  void* mSIRF_dataItems(const void* ptr_x);
EXPORTED_FUNCTION void* mSIRF_norm(const void* ptr_x);
EXPORTED_FUNCTION void*	mSIRF_dot(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION void* mSIRF_axpby(const PTR_FLOAT ptr_a, const void* ptr_x, const PTR_FLOAT ptr_b, const void* ptr_y);
EXPORTED_FUNCTION void* mSIRF_multiply(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION void* mSIRF_divide(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION void* mSIRF_write(const void* ptr, const char* filename);
EXPORTED_FUNCTION void* mSIRF_clone(void* ptr_x);
#ifndef CSIRF_FOR_MATLAB
}
#endif

#endif
