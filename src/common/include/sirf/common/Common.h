/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2026 Biomedical Research Foundation, Academy of Athens
Copyright 2026 University College London

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
\ingroup Common

\author Dimitra Kyriakopoulou
\author Kris Thielemans
*/

#ifndef SIRF_COMMON_H
#define SIRF_COMMON_H

#ifdef HAS_CUDA_RUNTIME_API
// cppcheck-suppress missingIncludeSystem
#include <cuda_runtime_api.h>
#endif

namespace sirf {

inline bool
pointer_supports_cuda_array_view(const void* ptr)
{
#ifdef HAS_CUDA_RUNTIME_API
    if (ptr == nullptr)
        return false;

    cudaPointerAttributes attrs{};
    const cudaError_t err = cudaPointerGetAttributes(&attrs, ptr);
    if (err != cudaSuccess)
        return false;
#if CUDART_VERSION >= 10000
    return attrs.type == cudaMemoryTypeManaged;
#else
    return attrs.memoryType == cudaMemoryTypeManaged && attrs.isManaged;
#endif
#else
    return false;
#endif
}

} // namespace sirf

#endif
