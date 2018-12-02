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

#include "data_handle.h"
#include "sirf/common/data_container.h"

using std::shared_ptr;
//#include "sirf/common/object_handle.inl"

using namespace sirf;

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
		s = z.real();
		return dataHandle(s);
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
		shared_ptr<DataContainer > sptr_z(x.new_data_container());
		sptr_z->axpby(ptr_a, x, ptr_b, y);
		return newObjectHandle<DataContainer >(sptr_z);
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
		shared_ptr<DataContainer > sptr_z(x.new_data_container());
		sptr_z->multiply(x, y);
		return newObjectHandle<DataContainer >(sptr_z);
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
		shared_ptr<DataContainer > sptr_z(x.new_data_container());
		sptr_z->divide(x, y);
		return newObjectHandle<DataContainer >(sptr_z);
	}
	CATCH;
}

