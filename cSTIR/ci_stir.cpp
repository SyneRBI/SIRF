#include "dh.h"
#include "stir.h"
#include "stir_p.h"
#include "stir_x.h"

#define CAST_PTR(T, X, Y) T* X = (T*)Y
#define NEW_SPTR(Base, X, Object) \
	boost::shared_ptr< Base >* X = new boost::shared_ptr< Base >(new Object)

static void*
unknownObject(const char* obj, const char* name, const char* file, int line)
{
	DataHandle* handle = new DataHandle;
	std::string error = "unknown ";
	error += obj;
	error += " '";
	error += name;
	error += "'";
	ExecutionStatus status(error.c_str(), file, line);
	handle->set(0, &status);
	return (void*)handle;
}

template<class Method>
void*
cSTIR_newReconstructionMethod(const char* parFile)
{
	try {
		if (strlen(parFile) > 0) {
			NEW_SPTR(Reconstruction<Image3DF>, ptr_sptr, Method(parFile));
			return newObjectHandle(ptr_sptr);
		}
		else {
			NEW_SPTR(Reconstruction<Image3DF>, ptr_sptr, Method);
			return newObjectHandle(ptr_sptr);
		}
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		DataHandle* handle = new DataHandle;
		handle->set(0, &status);
		return (void*)handle;
	}
	catch (...) {
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		DataHandle* handle = new DataHandle;
		handle->set(0, &status);
		return (void*)handle;
	}
}

extern "C"
void* cSTIR_newObject(const char* name)
{
	if (boost::iequals(name,
		"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
		return newObjectHandle
		< ObjectiveFunction3DF, PoissonLogLhLinModMeanProjData3DF >();
	if (boost::iequals(name, "ProjectorsUsingMatrix"))
		return newObjectHandle
		< ProjectorByBinPair, ProjectorPairUsingMatrix >();
	if (boost::iequals(name, "RayTracingMatrix"))
		return newObjectHandle
		< ProjMatrixByBin, RayTracingMatrix >();
	if (boost::iequals(name, "QuadraticPrior"))
		return newObjectHandle< Prior3DF, QuadPrior3DF >();
	if (boost::iequals(name, "TruncateToCylindricalFOVImageProcessor"))
		return newObjectHandle
		< DataProcessor3DF, CylindricFilter3DF >();
	if (boost::iequals(name, "EllipsoidalCylinder"))
		return newObjectHandle< Shape3D, EllipsoidalCylinder >();
	return unknownObject("object", name, __FILE__, __LINE__);
}

extern "C"
void* cSTIR_copyOfObject(void* ptr)
{
	CAST_PTR(anObjectHandle, ptr_obj, ptr);
	return (void*)ptr_obj->copy();
}

extern "C"
void cSTIR_deleteObject(void* ptr)
{
	if (!ptr)
		return;
	CAST_PTR(anObjectHandle, ptr_obj, ptr);
	delete ptr_obj;
}

extern "C"
void* cSTIR_setParameter
(void* ptr_s, const char* obj, const char* name, const void* ptr_v)
{
	CAST_PTR(DataHandle, hs, ptr_s);
	CAST_PTR(DataHandle, hv, ptr_v);
	if (boost::iequals(obj, "Shape"))
		return cSTIR_setShapeParameter(hs, name, hv);
	else if (boost::iequals(obj, "EllipsoidalCylinder"))
		return cSTIR_setEllipsoidalCylinderParameter(hs, name, hv);
	else if (boost::iequals(obj, "TruncateToCylindricalFOVImageProcessor"))
		return cSTIR_setTruncateToCylindricalFOVImageProcessorParameter
		(hs, name, hv);
	else if (boost::iequals(obj, "ProjectorsUsingMatrix"))
		return cSTIR_setProjectorsUsingMatrixParameter(hs, name, hv);
	else if (boost::iequals(obj, "RayTracingMatrix"))
		return cSTIR_setRayTracingMatrixParameter(hs, name, hv);
	else if (boost::iequals(obj, "GeneralisedPrior"))
		return cSTIR_setGeneralisedPriorParameter(hs, name, hv);
	else if (boost::iequals(obj, "QuadraticPrior"))
		return cSTIR_setQuadraticPriorParameter(hs, name, hv);
	else if (boost::iequals(obj, "GeneralisedObjectiveFunction"))
		return cSTIR_setGeneralisedObjectiveFunctionParameter(hs, name, hv);
	else if (boost::iequals(obj, "PoissonLogLikelihoodWithLinearModelForMean"))
		return cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanParameter
		(hs, name, hv);
	else if (boost::iequals(obj,
		"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
		return
		cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
		(hs, name, hv);
	else if (boost::iequals(obj, "Reconstruction"))
		return cSTIR_setReconstructionParameter(hs, name, hv);
	else if (boost::iequals(obj, "IterativeReconstruction"))
		return cSTIR_setIterativeReconstructionParameter(hs, name, hv);
	else if (boost::iequals(obj, "OSMAPOSL"))
		return cSTIR_setOSMAPOSLParameter(hs, name, hv);
	else if (boost::iequals(obj, "OSSPS"))
		return cSTIR_setOSSPSParameter(hs, name, hv);
	else
		return unknownObject("object", obj, __FILE__, __LINE__);
}

extern "C"
void* cSTIR_parameter(const void* ptr, const char* obj, const char* name) 
{
	CAST_PTR(DataHandle, handle, ptr);
	if (boost::iequals(obj, "Shape"))
		return cSTIR_shapeParameter(handle, name);
	if (boost::iequals(obj, "EllipsoidalCylinder"))
		return cSTIR_ellipsoidalCylinderParameter(handle, name);
	else if (boost::iequals(obj, "TruncateToCylindricalFOVImageProcessor"))
		return cSTIR_truncateToCylindricalFOVImageProcessorParameter
		(handle, name);
	if (boost::iequals(obj, "RayTracingMatrix"))
		return cSTIR_rayTracingMatrixParameter(handle, name);
	if (boost::iequals(obj, "ProjectorsUsingMatrix"))
		return cSTIR_projectorsUsingMatrixParameter(handle, name);
	if (boost::iequals(obj, "GeneralisedPrior"))
		return cSTIR_generalisedPriorParameter(handle, name);
	if (boost::iequals(obj, "GeneralisedObjectiveFunction"))
		return cSTIR_generalisedObjectiveFunctionParameter(handle, name);
	if (boost::iequals(obj,
		"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
		return
		cSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
		(handle, name);
	if (boost::iequals(obj, "IterativeReconstruction"))
		return cSTIR_iterativeReconstructionParameter(handle, name);
	if (boost::iequals(obj, "OSMAPOSL"))
		return cSTIR_OSMAPOSLParameter(handle, name);
	if (boost::iequals(obj, "OSSPS"))
		return cSTIR_OSSPSParameter(handle, name);
	return unknownObject("object", obj, __FILE__, __LINE__);
}

extern "C"
void* cSTIR_objectFromFile(const char* name, const char* filename)
{
	DataHandle* handle = new DataHandle;
	try {
		if (boost::iequals(name, "OSMAPOSLReconstruction"))
			return cSTIR_newReconstructionMethod
			<OSMAPOSLReconstruction<Image3DF> >
			(filename);
		if (boost::iequals(name, "OSSPSReconstruction"))
			return cSTIR_newReconstructionMethod
			<OSSPSReconstruction<Image3DF> >
			(filename);
		if (boost::iequals(name, "Image")) {
			sptrImage3DF* ptr_sptr = new sptrImage3DF
				(read_from_file<Image3DF>(filename));
			return newObjectHandle(ptr_sptr);
		}
		if (boost::iequals(name, "AcquisitionModelData")) {
			NEW(boost::shared_ptr<ProjData>, ptr_sptr);
			*ptr_sptr = ProjData::read_from_file(filename);
			return newObjectHandle(ptr_sptr);
		}
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	catch (...) {
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		handle->set(0, &status);
	}
	return (void*)handle;
}

extern "C"
void* cSTIR_setupObject(const char* obj, void* ptr_obj)
{
	CAST_PTR(DataHandle, ho, ptr_obj);
	DataHandle* handle = new DataHandle;
	try {
		bool status = 1;
		if (boost::iequals(obj, "GeneralisedPrior"))
			status = xSTIR_setupPrior(ho->data());
		else if (boost::iequals(obj, "GeneralisedObjectiveFunction"))
			status = xSTIR_setupObjectiveFunction(ho->data());
		if (status) {
			ExecutionStatus status("cSTIR_setupObject failed", __FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	catch (...) {
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		handle->set(0, &status);
	}
	return (void*)handle;
}

extern "C"
void* cSTIR_applyDataProcessor(const void* ptr_p, void* ptr_i)
{
	CAST_PTR(DataHandle, hp, ptr_p);
	CAST_PTR(DataHandle, hi, ptr_i);
	DataHandle* handle = new DataHandle;
	try {
		DataProcessor<Image3DF>& processor =
			objectFromHandle<DataProcessor<Image3DF> >(hp);
		Image3DF& image = objectFromHandle<Image3DF>(hi);
		processor.apply(image);
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	catch (...) {
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		handle->set(0, &status);
	}
	return (void*)handle;
}

extern "C"
void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i)
{
	CAST_PTR(DataHandle, hr, ptr_r);
	CAST_PTR(DataHandle, hi, ptr_i);
	DataHandle* handle = new DataHandle;
	try {
		sptrImage3DF& sptr_image = objectSptrFromHandle<Image3DF>(hi);
		Succeeded s = xSTIR_setupReconstruction(hr->data(), sptr_image);
		if (s != Succeeded::yes) {
			ExecutionStatus status("cSTIR_setupReconstruction failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	catch (...) {
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		handle->set(0, &status);
	}
	return (void*)handle;
}

extern "C"
void* cSTIR_runReconstruction(void* ptr_r, void* ptr_i) 
{
	CAST_PTR(DataHandle, hr, ptr_r);
	CAST_PTR(DataHandle, hi, ptr_i);
	DataHandle* handle = new DataHandle;
	try {
		Reconstruction<Image3DF>& recon =
			objectFromHandle< Reconstruction<Image3DF> >(hr);
		sptrImage3DF& sptr_image = objectSptrFromHandle<Image3DF>(hi);
		if (recon.reconstruct(sptr_image) != Succeeded::yes) {
			ExecutionStatus status("cSTIR_reconstruct failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	catch (...) {
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		handle->set(0, &status);
	}
	return (void*)handle;
}

extern "C"
void* cSTIR_updateReconstruction(void* ptr_r, void* ptr_i)
{
	CAST_PTR(DataHandle, hr, ptr_r);
	CAST_PTR(DataHandle, hi, ptr_i);
	DataHandle* handle = new DataHandle;
	try {
		Image3DF& image = objectFromHandle<Image3DF>(hi);
		xSTIR_updateReconstruction(hr->data(), image);
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	catch (...) {
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		handle->set(0, &status);
	}
	return (void*)handle;
}

extern "C"
void* cSTIR_voxels3DF
(int nx, int ny, int nz,
double sx, double sy, double sz,
double x, double y, double z)
{
	sptrVoxels3DF* sptr = new sptrVoxels3DF(
		new Voxels3DF(IndexRange3D(0, nz - 1,
			-(ny / 2), -(ny / 2) + ny - 1, -(nx / 2), -(nx / 2) + nx - 1),
			Coord3DF((float)z, (float)y, (float)x),
			Coord3DF((float)sz, (float)sy, (float)sx)));
	(*sptr)->fill(0);
	return newObjectHandle(sptr);
}

extern "C"
void* cSTIR_imageFromVoxels(void* ptr_v)
{
	CAST_PTR(DataHandle, hv, ptr_v);
	Voxels3DF& voxels = objectFromHandle<Voxels3DF>(hv);
	sptrImage3DF* sptr = new sptrImage3DF(voxels.clone());
	return newObjectHandle(sptr);
}

extern "C"
void* cSTIR_imageFromImage(void* ptr_i)
{
	CAST_PTR(DataHandle, hi, ptr_i);
	Image3DF& image = objectFromHandle<Image3DF>(hi);
	sptrImage3DF* sptr = new sptrImage3DF(image.clone());
	return newObjectHandle(sptr);
}

extern "C"
void* cSTIR_addShape(void* ptr_i, void* ptr_v, void* ptr_s, float v) 
{
	CAST_PTR(DataHandle, hi, ptr_i);
	CAST_PTR(DataHandle, hv, ptr_v);
	CAST_PTR(DataHandle, hs, ptr_s);

	Image3DF& image = objectFromHandle<Image3DF>(hi);
	Voxels3DF& voxels = objectFromHandle<Voxels3DF>(hv);
	Shape3D& shape = objectFromHandle<Shape3D>(hs);
	CartesianCoordinate3D<int> num_samples(1, 1, 1);
	voxels.fill(0);
	shape.construct_volume(voxels, num_samples);
	voxels *= v;
	image += voxels;

	return new DataHandle;
}

extern "C"
void cSTIR_fillImage(void* ptr_i, double v)
{
	CAST_PTR(DataHandle, hi, ptr_i);
	Image3DF& image = objectFromHandle<Image3DF>(hi);
	image.fill((float)v);
}

extern "C"
void cSTIR_getImageDimensions(const void* ptr, size_t pd) 
{
	CAST_PTR(DataHandle, handle, ptr);
	Image3DF& image = objectFromHandle<Image3DF>(handle);
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	int* dim = (int*)pd;
	if (!image.get_regular_range(min_indices, max_indices)) {
			//cout << "image is not regular" << endl;
		dim[0] = 0;
		dim[1] = 0;
		dim[2] = 0;
	}
	else {
		for (int i = 0; i < 3; i++)
			dim[i] = max_indices[i + 1] - min_indices[i + 1] + 1;
	}
}

extern "C"
void cSTIR_getImageData(const void* ptr, size_t pd) 
{
	CAST_PTR(DataHandle, handle, ptr);
	Image3DF& image = objectFromHandle<Image3DF>(handle);
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	double* data = (double*)pd;
	if (!image.get_regular_range(min_indices, max_indices)) {
			//cout << "image is not regular" << endl;
		return;
	}
	else {
		for (int z = min_indices[1], i = 0; z <= max_indices[1]; z++) {
			for (int y = min_indices[2]; y <= max_indices[2]; y++) {
				for (int x = min_indices[3]; x <= max_indices[3]; x++, i++) {
					data[i] = image[z][y][x];
					//data[i] = (*image)[z][y][x];
				}
			}
		}
	}
}

extern "C"
void* cSTIR_imagesDifference(void* first, void* second, int rimsize) 
{
	DataHandle* handle = new DataHandle;
	try {

		CAST_PTR(DataHandle, first_h, first);
		CAST_PTR(DataHandle, second_h, second);
		Image3DF& first_image = objectFromHandle<Image3DF>(first_h);
		Image3DF& second_image = objectFromHandle<Image3DF>(second_h);

		std::string explanation;
		if (!first_image.has_same_characteristics(second_image, explanation))
		{
			warning("input images do not have the same characteristics.\n%s",
				explanation.c_str());
			ExecutionStatus status(
				"input images do not have the same characteristics",
				__FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}

		if (rimsize >= 0)
		{
			truncate_rim(first_image, rimsize);
			truncate_rim(second_image, rimsize);
		}

		float reference_max = first_image.find_max();
		float reference_min = first_image.find_min();

		float amplitude = fabs(reference_max) > fabs(reference_min) ?
			fabs(reference_max) : fabs(reference_min);

		sptrImage3DF sptr(first_image.clone());
		Image3DF& image = *sptr.get();

		image -= second_image;
		const float max_error = image.find_max();
		const float min_error = image.find_min();

		float max_abs_error = fabs(min_error);
		if (max_error > max_abs_error)
			max_abs_error = max_error;

		double* result = (double*)malloc(sizeof(double));
		*result = max_abs_error / amplitude;
		handle->set(result, 0, GRAB);

	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	//catch (exception& e) {
	//	cout << "exception " << e.what() << endl;
	//}
	catch (...) {
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		handle->set(0, &status);
	}

	return (void*)handle;
}

