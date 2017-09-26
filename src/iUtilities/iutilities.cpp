/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

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

#include "shared_ptr.h"
#include "data_handle.h"

//using namespace SPTR_NAMESPACE;

extern "C" {

	void* newDataHandle() 
	{
		return (void*)new DataHandle;
	}
	void deleteDataHandle(void* ptr) 
	{
		if (ptr)
			delete (DataHandle*)ptr;
	}

	void deleteObject(void* ptr)
	{
		if (!ptr)
			return;
		CAST_PTR(DataHandle, ptr_obj, ptr);
		delete ptr_obj;
	}

	void* charDataHandle(const char* s) 
	{
		DataHandle* h = new DataHandle;
		size_t len = strlen(s);
		char* d = (char*)malloc(len + 1);
		//strcpy_s(d, len + 1, s);
		strcpy(d, s);
		h->set((void*)d, 0, GRAB);
		return (void*)h;
	}

	void* intDataHandle(int i)
	{
		return dataHandle<int>(i);
	}	
	void* floatDataHandle(float i) 
	{
		return dataHandle<float>(i);
	}
	void* doubleDataHandle(double i)
	{
		return dataHandle<double>(i);
	}

	char* charDataFromHandle(const void* ptr)
	{
		return charDataFromDataHandle((const DataHandle*)ptr);
	}
	int intDataFromHandle(const void* ptr)
	{
		return dataFromHandle<int>(ptr);
	}
	float floatDataFromHandle(const void* ptr) 
	{
		return dataFromHandle<float>(ptr);
	}
	float floatReDataFromHandle(const void* ptr)
	{
		std::complex<float> z = dataFromHandle<std::complex<float> >(ptr);
		return z.real();
	}
	float floatImDataFromHandle(const void* ptr)
	{
		std::complex<float> z = dataFromHandle<std::complex<float> >(ptr);
		return z.imag();
	}
	double doubleDataFromHandle(const void* ptr)
	{
		return dataFromHandle<double>(ptr);
	}
	double doubleReDataFromHandle(const void* ptr)
	{
		std::complex<double> z = dataFromHandle<std::complex<double> >(ptr);
		return z.real();
	}
	double doubleImDataFromHandle(const void* ptr)
	{
		std::complex<double> z = dataFromHandle<std::complex<double> >(ptr);
		return z.imag();
	}

	int executionStatus(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		return (ptr_h->status() ? 1 : 0);
	}

	const char* executionError(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		if (ptr_h->status())
			return ptr_h->status()->error();
		else
			return "";
	}

	const char* executionErrorFile(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		if (ptr_h->status())
			return ptr_h->status()->file();
		else
			return "";
	}

	int executionErrorLine(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		if (ptr_h->status())
			return ptr_h->status()->line();
		else
			return 0;
	}
}
