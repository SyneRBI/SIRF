#ifndef STIR_PARAMETERS_HANDLERS
#define STIR_PARAMETERS_HANDLERS

#include "dh.h"

#define NEW(T, X) T* X = new T
#define CAST_PTR(T, X, Y) T* X = (T*)Y
#define NEW_SPTR(T, X, Y) \
	boost::shared_ptr< T >* X = new boost::shared_ptr< T >(new Y)

template<class Base, class Object = Base>
Object&
objectFromHandle(const DataHandle* handle) {
	void* ptr = handle->data();
	if (ptr == 0)
		throw StirException("zero data pointer cannot be dereferenced", 
		__FILE__, __LINE__);
	CAST_PTR(boost::shared_ptr<Base>, sptr, ptr);
	if (is_null_ptr(*sptr))
		throw StirException("zero object pointer cannot be dereferenced",
		__FILE__, __LINE__);
	CAST_PTR(Object, object, sptr->get());
	return *object;
}

template<class T>
void*
sptrDataHandle(boost::shared_ptr<T> sptr) {
	NEW(boost::shared_ptr<T>, ptr_sptr);
	*ptr_sptr = sptr;
	NEW(DataHandle, handle);
	handle->set((void*)ptr_sptr, 0);
	return (void*)handle;
}

template<class T>
boost::shared_ptr<T>
sptrDataFromHandle(const DataHandle* handle) {
	return *(boost::shared_ptr<T>*)handle->data();
}

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