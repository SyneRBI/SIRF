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

#include <stddef.h>

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


// MR Simulation
void* cDS_MRDynamicSimulation(const void* ptr_labels, const char* fname_xml);
void* cDS_simulateData(void* ptr_sim);
void* cDS_writeSimulationResults(const void* ptr_sim, const char* fname_with_ext);
void* cDS_saveMotionGroundTruth(const void* ptr_sim);

void* cDS_setAcquisitionTemplateData(void* ptr_sim, const void* ptr_acqs);
void* cDS_setContrastTemplateData(void* ptr_sim, const void* ptr_acqs);
void* cDS_setCoilmaps(void* ptr_sim, const void* ptr_csm);

void* cDS_setSNR(void* ptr_sim, float const SNR);
void* cDS_setNoiseLabel(void* ptr_sim, int const label);

void* cDS_setOffsetTransformation(void* ptr_sim, const void* ptr_trafo);

void* cDS_addMRMotionDynamic(void* ptr_sim, void* ptr_dyn); 
void* cDS_addMRContrastDynamic(void* ptr_sim, void* ptr_dyn); 
void* cDS_addExternalContrastDynamic(void* ptr_sim, void* ptr_dyn);
// Dynamics
void* cDS_setDynamicSignal(void* ptr_dyn, PTR_FLOAT ptr_time, PTR_FLOAT ptr_signal, int const num_points);
void* cDS_addMRDisplacementField(void* ptr_dyn, const void* ptr_dvf);
void* cDS_setMRGroundTruthFolderName(void* ptr_dyn, const char* fpath_output_prefix);
void* cDS_setCyclicality(void* ptr_dyn, bool const cyc);


void* cDS_setMRAcquisitions(void* ptr_dyn, void* ptr_ad);

// MR Dynamics
void* cDS_MRMotionDynamic( int const num_states );

void* cDS_MRContrastDynamic( int const num_states );
void* cDS_addDynamicLabel(void* ptr_contrast_dyn, int const label);
void* cDS_setMRParameterExtremes(void* ptr_contrast_dyn, void* ptr_tissueparameter_0, void* ptr_tissueparameter_1);

void* cDS_ExternalMRContrastDynamic( void );
void* cDS_appendExternalTissueSignal(void* ptr_dyn, int const num_points, PTR_INT ptr_labels, PTR_FLOAT ptr_sig);


// Tissue Parameters
void* cDS_getTissueParameter(const void* ptr_sim, const unsigned int label);
void* cDS_setT1Value(void* ptr_tissue_parameter, float const T1_ms);


#ifndef CSIMULATION_FOR_MATLAB
}
#endif

#endif
