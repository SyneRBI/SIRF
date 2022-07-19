/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC
Copyright 2017 - 2020 University College London

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

#ifndef STIR_PARAMETERS_HANDLERS
#define STIR_PARAMETERS_HANDLERS

#include "sirf/iUtilities/DataHandle.h"

namespace sirf {

	void*
		cSTIR_setImageDataParameter(void* hp, const char* name, const void* hv);

	void*
		cSTIR_ImageDataParameter(void* hp, const char* name);

	void*
		cSTIR_setListmodeToSinogramsParameter(void* hp, const char* name, const void* hv);

	void*
		cSTIR_setShapeParameter(void* hp, const char* name, const void* hv);

	void*
		cSTIR_shapeParameter(const DataHandle* handle, const char* name);

	void*
		cSTIR_setEllipsoidalCylinderParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_ellipsoidalCylinderParameter(const DataHandle* handle, const char* name);

	void*
		cSTIR_setRayTracingMatrixParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_rayTracingMatrixParameter(const DataHandle* handle, const char* name);

	void*
		cSTIR_setSPECTUBMatrixParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_SPECTUBMatrixParameter(const DataHandle* handle, const char* name);

	void*
		cSTIR_setAcquisitionModelParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_AcquisitionModelParameter(DataHandle* hm, const char* name);

	void*
		cSTIR_setAcqModUsingMatrixParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

#ifdef STIR_WITH_NiftyPET_PROJECTOR
    void*
        cSTIR_setAcqModUsingNiftyPETParameter
        (DataHandle* hp, const char* name, const DataHandle* hv);
#endif

	void*
		cSTIR_acqModUsingMatrixParameter(DataHandle* hp, const char* name);

	void*
		cSTIR_setTruncateToCylindricalFOVImageProcessorParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_truncateToCylindricalFOVImageProcessorParameter
		(const DataHandle* handle, const char* name);

	void*
		cSTIR_setSeparableGaussianImageFilterParameter
		(void* hp, const char* name, const void* hv);

	void*
		cSTIR_setGeneralisedPriorParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_generalisedPriorParameter(const DataHandle* handle, const char* name);

	void*
		cSTIR_setQuadraticPriorParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);
	void*
		cSTIR_QuadraticPriorParameter
		(DataHandle* hp, const char* name);

	void*
		cSTIR_setRelativeDifferencePriorParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);
	void*
		cSTIR_RelativeDifferencePriorParameter
		(DataHandle* hp, const char* name);

	void*
		cSTIR_setPLSPriorParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_PLSPriorParameter
		(DataHandle* hp, const char* name);

	void*
		cSTIR_generalisedObjectiveFunctionParameter
		(const DataHandle* handle, const char* name);

        void*
                cSTIR_setScatterSimulatorParameter
                (const DataHandle* hp, const char* name, const DataHandle* hv);
        void*
                cSTIR_setScatterEstimatorParameter
                (const DataHandle *hp, const char* name, const DataHandle* hv);
        void*
                cSTIR_ScatterEstimatorParameter
                (DataHandle* hp, const char* name);

	void*
		cSTIR_setGeneralisedObjectiveFunctionParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
		(const DataHandle* handle, const char* name);

	void*
		cSTIR_setReconstructionParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_setIterativeReconstructionParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_iterativeReconstructionParameter(const DataHandle* handle, const char* name);

	void*
		cSTIR_setOSMAPOSLParameter
		(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_OSMAPOSLParameter(const DataHandle* handle, const char* name);
#ifdef USE_HKEM
	void*
		cSTIR_setKOSMAPOSLParameter(DataHandle* hp, const char* name, const DataHandle* hv);
#endif
	void*
		cSTIR_setOSSPSParameter(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_OSSPSParameter(const DataHandle* handle, const char* name);

	void*
		cSTIR_setFBP2DParameter(DataHandle* hp, const char* name, const DataHandle* hv);

	void*
		cSTIR_FBP2DParameter(DataHandle* hp, const char* name);

}

#endif
