#ifndef STIR_PARAMETERS_HANDLERS
#define STIR_PARAMETERS_HANDLERS

#include "data_handle.h"

void*
cSTIR_setShapeParameter(DataHandle* hp, const char* name, const DataHandle* hv);

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
cSTIR_setProjectorsUsingMatrixParameter
(DataHandle* hp, const char* name, const DataHandle* hv);

void*
cSTIR_projectorsUsingMatrixParameter
(const DataHandle* handle, const char* name);

void*
cSTIR_setTruncateToCylindricalFOVImageProcessorParameter
(DataHandle* hp, const char* name, const DataHandle* hv);

void*
cSTIR_truncateToCylindricalFOVImageProcessorParameter
(const DataHandle* handle, const char* name);

void*
cSTIR_setGeneralisedPriorParameter
(DataHandle* hp, const char* name, const DataHandle* hv);

void*
cSTIR_generalisedPriorParameter(const DataHandle* handle, const char* name);

void*
cSTIR_setQuadraticPriorParameter
(DataHandle* hp, const char* name, const DataHandle* hv);

void*
cSTIR_generalisedObjectiveFunctionParameter
(const DataHandle* handle, const char* name);

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

void*
cSTIR_setOSSPSParameter(DataHandle* hp, const char* name, const DataHandle* hv);

void*
cSTIR_OSSPSParameter(const DataHandle* handle, const char* name);

#endif