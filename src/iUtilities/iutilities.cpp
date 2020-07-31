/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2018 Rutherford Appleton Laboratory STFC
Copyright 2019 University College London

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

/*!
\file
\ingroup C Interface to C++ Objects
\brief C interface to DataHandle.

Defines C functions handling DataHandle objects.
\author Evgueni Ovtchinnikov
\author SyneRBI
*/

#include <stdint.h>
#include <complex>

#include "sirf/iUtilities/DataHandle.h"

extern "C" {

	void* newDataHandle() // C constructor
	{
		return (void*)new DataHandle;
	}
	void deleteDataHandle(void* ptr) // C destructor
	{
		if (ptr)
			delete (DataHandle*)ptr;
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
    bool boolDataFromHandle(const void* ptr)
    {
        return dataFromHandle<bool>(ptr);
    }
	int intDataFromHandle(const void* ptr)
	{
		return dataFromHandle<int>(ptr);
	}
	int intDataItemFromHandle(const void* ptr, int i)
	{
		int* arr = dataFromHandle<int*>(ptr);
		return arr[i];
	}
	int uint16DataItemFromHandle(const void* ptr, int i)
	{
		uint16_t* arr = dataFromHandle<uint16_t*>(ptr);
		return arr[i];
	}
	int uint32DataItemFromHandle(const void* ptr, int i)
	{
		uint32_t* arr = dataFromHandle<uint32_t*>(ptr);
		return arr[i];
	}
	int uint64DataItemFromHandle(const void* ptr, int i)
	{
		uint64_t* arr = dataFromHandle<uint64_t*>(ptr);
		return arr[i];
	}
	float floatDataFromHandle(const void* ptr)
	{
		return dataFromHandle<float>(ptr);
	}
	float floatDataItemFromHandle(const void* ptr, int i)
	{
		float* arr = dataFromHandle<float*>(ptr);
		return arr[i];
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
			return ptr_h->status()->error().c_str();
		else
			return "";
	}

	const char* executionErrorFile(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		if (ptr_h->status())
			return ptr_h->status()->file().c_str();
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
