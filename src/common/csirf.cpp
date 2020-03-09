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

#include <complex>
#include <memory>

#include "sirf/iUtilities/DataHandle.h"
#include "sirf/common/DataContainer.h"
#include "sirf/common/ImageData.h"

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
		DataContainer& x =
			objectFromHandle<DataContainer >(h_x);
		return dataHandle(x.items());
	}
	CATCH;
}

extern "C"
void*
cSIRF_norm(const void* ptr_x)
{
	try {
		DataContainer& x =
			objectFromHandle<DataContainer >(ptr_x);
		return dataHandle(x.norm());
	}
	CATCH;
}

extern "C"
void*
cSIRF_dot(const void* ptr_x, const void* ptr_y)
{
	try {
		DataContainer& x =
			objectFromHandle<DataContainer >(ptr_x);
		DataContainer& y =
			objectFromHandle<DataContainer >(ptr_y);
		float s;
		std::complex<float> z(0.0, 0.0);
		x.dot(y, &z);
		//s = z.real();
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
		DataContainer& x =
			objectFromHandle<DataContainer >(ptr_x);
		DataContainer& y =
			objectFromHandle<DataContainer >(ptr_y);
		void* h = x.new_data_container_handle();
		DataContainer& z = objectFromHandle<DataContainer>(h);
		z.axpby(ptr_a, x, ptr_b, y);
		return h;
		//shared_ptr<DataContainer > sptr_z(x.new_data_container());
		//sptr_z->axpby(ptr_a, x, ptr_b, y);
		//return newObjectHandle<DataContainer >(sptr_z);
	}
	CATCH;
}

extern "C"
void*
cSIRF_multiply(const void* ptr_x, const void* ptr_y)
{
	try {
		DataContainer& x =
			objectFromHandle<DataContainer >(ptr_x);
		DataContainer& y =
			objectFromHandle<DataContainer >(ptr_y);
		void* h = x.new_data_container_handle();
		DataContainer& z = objectFromHandle<DataContainer>(h);
		z.multiply(x, y);
		return h;
		//shared_ptr<DataContainer > sptr_z(x.new_data_container());
		//sptr_z->multiply(x, y);
		//return newObjectHandle<DataContainer >(sptr_z);
	}
	CATCH;
}

extern "C"
void*
cSIRF_divide(const void* ptr_x, const void* ptr_y)
{
	try {
		DataContainer& x =
			objectFromHandle<DataContainer >(ptr_x);
		DataContainer& y =
			objectFromHandle<DataContainer >(ptr_y);
		void* h = x.new_data_container_handle();
		DataContainer& z = objectFromHandle<DataContainer>(h);
		z.divide(x, y);
		return h;
		//shared_ptr<DataContainer > sptr_z(x.new_data_container());
		//sptr_z->divide(x, y);
		//return newObjectHandle<DataContainer >(sptr_z);
	}
	CATCH;
}

extern "C"
void*
cSIRF_write(const void* ptr, const char* filename)
{
	try {
		DataContainer& data =
			objectFromHandle<DataContainer >(ptr);
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
		DataContainer& x =
			objectFromHandle<DataContainer >(ptr_x);
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
	ImageData& id = objectFromHandle<ImageData>(ptr_im);
	ImageData& id_src = objectFromHandle<ImageData>(ptr_src);
	id.fill(id_src);
	return new DataHandle;
}

extern "C"
void* cSIRF_ImageData_reorient(void* im_ptr, void *geom_info_ptr)
{
    try {
        ImageData& id = objectFromHandle<ImageData>(im_ptr);
        VoxelisedGeometricalInfo3D geom_info =
                objectFromHandle<VoxelisedGeometricalInfo3D>(geom_info_ptr);
        id.reorient(geom_info);
        return new DataHandle;
    }
    CATCH
}

extern "C"
void*
cSIRF_ImageData_get_geom_info(const void* ptr_im)
{
    const ImageData& id = objectFromHandle<const ImageData>(ptr_im);
    return newObjectHandle(id.get_geom_info_sptr());
}

extern "C"
void*
cSIRF_GeomInfo_print(const void* ptr_geom)
{
    const VoxelisedGeometricalInfo3D &geom_info =
            objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
    geom_info.print_info();
    return new DataHandle;
}

extern "C"
void*
cSIRF_GeomInfo_get_offset(const void* ptr_geom, void* ptr_arr)
{
    const VoxelisedGeometricalInfo3D &geom_info =
            objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
    const VoxelisedGeometricalInfo3D::Offset offset =
            geom_info.get_offset();
    float *data = (float*)ptr_arr;
    for (unsigned i=0; i<3; ++i)
        data[i] = offset[i];
    return new DataHandle;
}

extern "C"
void*
cSIRF_GeomInfo_get_spacing(const void* ptr_geom, void* ptr_arr)
{
    const VoxelisedGeometricalInfo3D &geom_info =
            objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
    const VoxelisedGeometricalInfo3D::Spacing spacing =
            geom_info.get_spacing();
    float *data = (float*)ptr_arr;
    for (unsigned i=0; i<3; ++i)
        data[i] = spacing[i];
    return new DataHandle;
}

extern "C"
void*
cSIRF_GeomInfo_get_size(const void* ptr_geom, void* ptr_arr)
{
    const VoxelisedGeometricalInfo3D &geom_info =
            objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
    const VoxelisedGeometricalInfo3D::Size size =
            geom_info.get_size();
    int *data = (int*)ptr_arr;
    for (unsigned i=0; i<3; ++i)
        data[i] = size[i];
    return new DataHandle;
}

extern "C"
void*
cSIRF_GeomInfo_get_direction_matrix(const void* ptr_geom, void* ptr_arr)
{
    const VoxelisedGeometricalInfo3D &geom_info =
            objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
    const VoxelisedGeometricalInfo3D::DirectionMatrix dm =
            geom_info.get_direction();
    float *data = (float*)ptr_arr;
    for (unsigned i=0; i<3; ++i)
        for (unsigned j=0; j<3; ++j)
        data[i*3+j] = dm[i][j];
    return new DataHandle;
}

extern "C"
void*
cSIRF_GeomInfo_get_index_to_physical_point_matrix(const void* ptr_geom, void* ptr_arr)
{
    const VoxelisedGeometricalInfo3D &geom_info =
            objectFromHandle<const VoxelisedGeometricalInfo3D>(ptr_geom);
    const VoxelisedGeometricalInfo3D::TransformMatrix tm =
            geom_info.calculate_index_to_physical_point_matrix();
    float *data = (float*)ptr_arr;
    for (unsigned i=0; i<4; ++i)
        for (unsigned j=0; j<4; ++j)
        data[i*3+j] = tm[i][j];
    return new DataHandle;
}