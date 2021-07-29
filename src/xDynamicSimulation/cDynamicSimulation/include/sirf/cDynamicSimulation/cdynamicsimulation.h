/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2021 Physikalisch-Technische Bundesanstalt (PTB)

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

#ifndef SIMULATION_C_INTERFACE
#define SIMULATION_C_INTERFACE

#ifndef CSIMULATION_FOR_MATLAB
#define PTR_INT size_t
#define PTR_FLOAT size_t
#define PTR_DOUBLE size_t
extern "C" {
#else
#define PTR_INT int*
#define PTR_FLOAT float*
#define PTR_DOUBLE double*
#endif

#include "sirf/cDynamicSimulation/cdynamicsimulation.h"

// dynamic simulation 
    void* cDS_MRDynamicSimulation(void* ptr_labels, const char* fname_xml);

#ifndef CSIMULATION_FOR_MATLAB
}
#endif

#endif
