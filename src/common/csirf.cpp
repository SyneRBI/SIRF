/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2020 University College London.

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

#include <complex>
#include <memory>

#include "sirf/iUtilities/DataHandle.h"
#include "sirf/common/DataContainer.h"
#include "sirf/common/iequals.h"
#include "sirf/common/ImageData.h"
#include "sirf/Syn/utilities.h"
#include "sirf/common/deprecate.h"

using namespace sirf;

#define NEW_OBJECT_HANDLE(T) new ObjectHandle<T >(shared_ptr<T >(new T))
#define SPTR_FROM_HANDLE(Object, X, H) \
	shared_ptr<Object> X; getObjectSptrFromHandle<Object>(H, X);


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

//default constructors
extern "C"
void* cSIRF_newObject(const char* name)
{
	try {
        if (strcmp(name, "DataHandleVector") == 0)
            return newObjectHandle(std::shared_ptr<DataHandleVector>(new DataHandleVector));
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}


extern "C"
void*
cSIRF_dataItems(const void* ptr_x)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		auto const& x =	objectFromHandle<DataContainer>(h_x);
		return dataHandle(x.items());
	}
	CATCH;
}

extern "C"
void*
cSIRF_isComplex(const void* ptr_x)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		auto const& x =	objectFromHandle<DataContainer>(h_x);
		return dataHandle<int>(x.is_complex());
	}
	CATCH;
}

extern "C"
void*
cSIRF_bits(const void* ptr_x)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		auto const& x = objectFromHandle<DataContainer>(h_x);
		return dataHandle<int>(x.bits());
	}
	CATCH;
}

extern "C"
void*
cSIRF_conjugate(void* ptr)
{
	try {
		auto& x = objectFromHandle<DataContainer>(ptr);
		x.conjugate();
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_conjugated(void* ptr_x)
{
	try {
		auto const& x =	objectFromHandle<DataContainer>(ptr_x);
		std::shared_ptr<DataContainer> sptr(x.clone());
		sptr->conjugate();
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void*
cSIRF_norm(const void* ptr_x)
{
	try {
		auto const& x =	objectFromHandle<DataContainer>(ptr_x);
		return dataHandle(x.norm());
	}
	CATCH;
}

extern "C"
void*
cSIRF_dot(const void* ptr_x, const void* ptr_y)
{
	try {
		auto const& x =	objectFromHandle<DataContainer>(ptr_x);
		auto const& y =	objectFromHandle<DataContainer>(ptr_y);
		std::complex<float> z(0.0, 0.0);
		x.dot(y, &z);
		return dataHandle(z);
	}
	CATCH;
}

extern "C"
void*
cSIRF_integral(const void* ptr_x)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		std::complex<float> z(0.0, 0.0);
		x.sum(&z);
		return dataHandle(z);
	}
	CATCH;
}

extern "C"
void*
cSIRF_max(const void* ptr_x)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		std::complex<float> z(0.0, 0.0);
		x.max(&z);
		return dataHandle(z);
	}
	CATCH;
}

extern "C"
void*
cSIRF_axpby(
	const void* ptr_a, const void* ptr_x,
	const void* ptr_b, const void* ptr_y
) {
	try {
		auto const& x =	objectFromHandle<DataContainer>(ptr_x);
		auto const& y =	objectFromHandle<DataContainer>(ptr_y);
		void* h = x.new_data_container_handle();
		auto& z = objectFromHandle<DataContainer>(h);
		z.xapyb(x, ptr_a, y, ptr_b);
		return h;
	}
	CATCH;
}

extern "C"
void*
cSIRF_axpbyAlt(
	const void* ptr_a, const void* ptr_x,
	const void* ptr_b, const void* ptr_y,
	void* ptr_z
) {
	try {
		auto const& x =	objectFromHandle<DataContainer>(ptr_x);
		auto const& y =	objectFromHandle<DataContainer>(ptr_y);
		auto& z = objectFromHandle<DataContainer>(ptr_z);
		z.xapyb(x, ptr_a, y, ptr_b);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_xapyb(
	const void* ptr_x, const void* ptr_a,
	const void* ptr_y, const void* ptr_b
) {
	try {
		auto const& x =	objectFromHandle<DataContainer>(ptr_x);
		auto const& a =	objectFromHandle<DataContainer>(ptr_a);
		auto const& y =	objectFromHandle<DataContainer>(ptr_y);
		auto const& b =	objectFromHandle<DataContainer>(ptr_b);
		void* h = x.new_data_container_handle();
		auto& z = objectFromHandle<DataContainer>(h);
		z.xapyb(x, a, y, b);
		return h;
	}
	CATCH;
}

extern "C"
void*
cSIRF_xapybAlt(
	const void* ptr_x, const void* ptr_a,
	const void* ptr_y, const void* ptr_b,
	void* ptr_z
) {
	try {
		auto const& x =	objectFromHandle<DataContainer>(ptr_x);
		auto const& a =	objectFromHandle<DataContainer>(ptr_a);
		auto const& y =	objectFromHandle<DataContainer>(ptr_y);
		auto const& b =	objectFromHandle<DataContainer>(ptr_b);
		auto& z = objectFromHandle<DataContainer>(ptr_z);
		z.xapyb(x, a, y, b);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_XapYB(
	const void* ptr_x, const void* ptr_a,
	const void* ptr_y, const void* ptr_b
) {
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		auto const& y = objectFromHandle<DataContainer>(ptr_y);
		auto const& b = objectFromHandle<DataContainer>(ptr_b);
		void* h = x.new_data_container_handle();
		auto& z = objectFromHandle<DataContainer>(h);
		z.xapyb(x, ptr_a, y, b);
		return h;
	}
	CATCH;
}

extern "C"
void*
cSIRF_XapYBAlt(
	const void* ptr_x, const void* ptr_a,
	const void* ptr_y, const void* ptr_b,
	void* ptr_z
) {
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		auto const& y = objectFromHandle<DataContainer>(ptr_y);
		auto const& b = objectFromHandle<DataContainer>(ptr_b);
		auto& z = objectFromHandle<DataContainer>(ptr_z);
		z.xapyb(x, ptr_a, y, b);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_add(const void* ptr_x, const void* ptr_y, const void* ptr_z)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		auto& z = objectFromHandle<DataContainer>(ptr_z);
		z.add(x, ptr_y);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_sum(const void* ptr_x, const void* ptr_y)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		void* h = x.new_data_container_handle();
		auto& z = objectFromHandle<DataContainer>(h);
		z.add(x, ptr_y);
		return h;
	}
	CATCH;
}

extern "C"
void*
cSIRF_binary(const void* ptr_x, const void* ptr_y, const char* f)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		auto const& y = objectFromHandle<DataContainer>(ptr_y);
		void* h = x.new_data_container_handle();
		auto& z = objectFromHandle<DataContainer>(h);
		if (sirf::iequals(f, "power"))
			z.power(x, y);
		else if (sirf::iequals(f, "multiply"))
			z.multiply(x, y);
		else if (sirf::iequals(f, "divide"))
			z.divide(x, y);
		else if (sirf::iequals(f, "maximum"))
			z.maximum(x, y);
		else if (sirf::iequals(f, "minimum"))
			z.minimum(x, y);
		else
			return unknownObject("function", f, __FILE__, __LINE__);
		return h;
	}
	CATCH;
}

extern "C"
void*
cSIRF_compute_binary(const void* ptr_x, const void* ptr_y, const char* f, const void* ptr_z)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		auto& z = objectFromHandle<DataContainer>(ptr_z);
		auto const& y = objectFromHandle<DataContainer>(ptr_y);
		if (sirf::iequals(f, "power"))
			z.power(x, y);
		else if (sirf::iequals(f, "multiply"))
			z.multiply(x, y);
		else if (sirf::iequals(f, "divide"))
			z.divide(x, y);
		else if (sirf::iequals(f, "maximum"))
			z.maximum(x, y);
		else if (sirf::iequals(f, "minimum"))
			z.minimum(x, y);
		else
			return unknownObject("function", f, __FILE__, __LINE__);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_semibinary(const void* ptr_x, const void* ptr_y, const char* f)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		void* h = x.new_data_container_handle();
		auto& z = objectFromHandle<DataContainer>(h);
		if (sirf::iequals(f, "power"))
			z.power(x, ptr_y);
		else if (sirf::iequals(f, "multiply"))
			z.multiply(x, ptr_y);
		else if (sirf::iequals(f, "maximum"))
			z.maximum(x, ptr_y);
		else if (sirf::iequals(f, "minimum"))
			z.minimum(x, ptr_y);
		else
			return unknownObject("function", f, __FILE__, __LINE__);
		return h;
	}
	CATCH;
}

extern "C"
void*
cSIRF_compute_semibinary(const void* ptr_x, const void* ptr_y, const char* f, const void* ptr_z)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		auto& z = objectFromHandle<DataContainer>(ptr_z);
		if (sirf::iequals(f, "power"))
			z.power(x, ptr_y);
		else if (sirf::iequals(f, "multiply"))
			z.multiply(x, ptr_y);
		else if (sirf::iequals(f, "maximum"))
			z.maximum(x, ptr_y);
		else if (sirf::iequals(f, "minimum"))
			z.minimum(x, ptr_y);
		else
			return unknownObject("function", f, __FILE__, __LINE__);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_unary(const void* ptr_x, const char* f)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		void* h = x.new_data_container_handle();
		auto& z = objectFromHandle<DataContainer>(h);
		if (sirf::iequals(f, "exp"))
			z.exp(x);
		else if (sirf::iequals(f, "log"))
			z.log(x);
		else if (sirf::iequals(f, "sqrt"))
			z.sqrt(x);
		else if (sirf::iequals(f, "sign"))
			z.sign(x);
		else if (sirf::iequals(f, "abs"))
			z.abs(x);
		else
			return unknownObject("function", f, __FILE__, __LINE__);
		return h;
	}
	CATCH;
}

extern "C"
void*
cSIRF_compute_unary(const void* ptr_x, const char* f, const void* ptr_z)
{
	try {
		auto const& x = objectFromHandle<DataContainer>(ptr_x);
		auto& z = objectFromHandle<DataContainer>(ptr_z);
		if (sirf::iequals(f, "exp"))
			z.exp(x);
		else if (sirf::iequals(f, "log"))
			z.log(x);
		else if (sirf::iequals(f, "sqrt"))
			z.sqrt(x);
		else if (sirf::iequals(f, "sign"))
			z.sign(x);
		else if (sirf::iequals(f, "abs"))
			z.abs(x);
		else
			return unknownObject("function", f, __FILE__, __LINE__);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_write(const void* ptr, const char* filename)
{
	try {
		auto const& data = objectFromHandle<DataContainer>(ptr);
		data.write(filename);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* 
cSIRF_clone(void* ptr_x)
{
	try {
		auto const& x =	objectFromHandle<DataContainer>(ptr_x);
		std::shared_ptr<DataContainer> sptr(x.clone());
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void*
cSIRF_DataHandleVector_push_back(void* self, void* to_append)
{
    DataHandleVector& vec = objectFromHandle<DataHandleVector>(self);
    vec.push_back(to_append);
    return new DataHandle;
}

extern "C"
void*
cSIRF_fillImageFromImage(void* ptr_im, const void* ptr_src)
{
	try {
		auto& id = objectFromHandle<ImageData>(ptr_im);
		auto const& id_src = objectFromHandle<ImageData>(ptr_src);
		id.fill(id_src);
		return new DataHandle;
	}
    CATCH;
}

extern "C"
void*
cSIRF_readImageData(const char* file, const char* eng, int verb)
{
	try {
		ImageDataWrap idw(file, eng, verb);
		std::shared_ptr<ImageData> sptr_id = idw.data_sptr();
		return newObjectHandle<ImageData>(sptr_id);
	}
	CATCH;
}

extern "C"
void* 
cSIRF_equalImages(const void* ptr_im_a, const void* ptr_im_b)
{
	try {
		auto const& id_a = objectFromHandle<ImageData>(ptr_im_a);
		auto const& id_b = objectFromHandle<ImageData>(ptr_im_b);
		int same = (id_a == id_b);
		return dataHandle(same);
	}
    CATCH;
}

extern "C"
void* 
cSIRF_ImageData_reorient(void* im_ptr, void *geom_info_ptr)
{
    try {
        auto& id = objectFromHandle<ImageData>(im_ptr);
        VoxelisedGeometricalInfo3D geom_info =
                objectFromHandle<VoxelisedGeometricalInfo3D>(geom_info_ptr);
        id.reorient(geom_info);
        return new DataHandle;
    }
    CATCH;
}

extern "C"
void*
cSIRF_ImageData_get_geom_info(const void* ptr_im)
{
	try {
		const auto& id = objectFromHandle<const ImageData>(ptr_im);
		return newObjectHandle(id.get_geom_info_sptr());
	}
	CATCH;
}

extern "C"
void*
cSIRF_GeomInfo_get(const void* ptr_geom)
{
	try {
		const auto& geom_info =
			objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
		return charDataHandleFromCharData(geom_info.get_info().c_str());
	}
	CATCH;
}

extern "C"
void*
cSIRF_GeomInfo_get_offset(const void* ptr_geom, void* ptr_arr)
{
	try {
		const auto& geom_info =
			objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
		const VoxelisedGeometricalInfo3D::Offset offset =
			geom_info.get_offset();
		float *data = (float*)ptr_arr;
		for (unsigned i = 0; i < 3; ++i)
			data[i] = offset[i];
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_GeomInfo_get_spacing(const void* ptr_geom, void* ptr_arr)
{
	try {
		const auto& geom_info =
			objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
		const VoxelisedGeometricalInfo3D::Spacing spacing =
			geom_info.get_spacing();
		float *data = (float*)ptr_arr;
		for (unsigned i = 0; i < 3; ++i)
			data[i] = spacing[i];
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_GeomInfo_get_size(const void* ptr_geom, void* ptr_arr)
{
	try {
		const auto& geom_info =
			objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
		const VoxelisedGeometricalInfo3D::Size size =
			geom_info.get_size();
		int *data = (int*)ptr_arr;
		for (unsigned i = 0; i < 3; ++i)
			data[i] = size[i];
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_GeomInfo_get_direction_matrix(const void* ptr_geom, void* ptr_arr)
{
	try {
		const auto& geom_info =
			objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
		const VoxelisedGeometricalInfo3D::DirectionMatrix dm =
			geom_info.get_direction();
		float *data = (float*)ptr_arr;
		for (unsigned i = 0; i < 3; ++i)
			for (unsigned j = 0; j < 3; ++j)
				data[i * 3 + j] = dm[i][j];
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSIRF_GeomInfo_get_index_to_physical_point_matrix(const void* ptr_geom, void* ptr_arr)
{
	try {
		const auto& geom_info =
			objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
		const VoxelisedGeometricalInfo3D::TransformMatrix tm =
			geom_info.calculate_index_to_physical_point_matrix();
		float *data = (float*)ptr_arr;
		for (unsigned i = 0; i < 4; ++i)
			for (unsigned j = 0; j < 4; ++j)
				data[i * 4 + j] = tm[i][j];
		return new DataHandle;
	}
	CATCH;
}
