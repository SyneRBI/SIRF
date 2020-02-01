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

#ifndef INTERFACE_UTILITIES
#define INTERFACE_UTILITIES

// Deprecation function. With C++14, could use [[deprecated("some message")]]
#if defined(__GNUC__) || defined(__clang__)
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

#ifndef IUTILITIES_FOR_MATLAB
extern "C" {
#endif
	void* newDataHandle();
	void deleteDataHandle(void* ptr);
	void* charDataHandle(const char* s);
	void* intDataHandle(int i);
	void* floatDataHandle(float i);
	void* doubleDataHandle(double i);
	char* charDataFromHandle(const void* ptr);
    bool boolDataFromHandle(const void* ptr);
	int intDataFromHandle(const void* ptr);
	int intDataItemFromHandle(const void* ptr, int i);
	int uint16DataItemFromHandle(const void* ptr, int i);
	int uint32DataItemFromHandle(const void* ptr, int i);
	int uint64DataItemFromHandle(const void* ptr, int i);
	float floatDataFromHandle(const void* ptr);
	float floatDataItemFromHandle(const void* ptr, int i);
	float floatReDataFromHandle(const void* ptr);
	float floatImDataFromHandle(const void* ptr);
	double doubleDataFromHandle(const void* ptr);
	double doubleReDataFromHandle(const void* ptr);
	double doubleImDataFromHandle(const void* ptr);
	int executionStatus(const void* ptr);
	const char* executionError(const void* ptr);
	const char* executionErrorFile(const void* ptr);
	int executionErrorLine(const void* ptr);
#ifndef IUTILITIES_FOR_MATLAB
}
#endif

#endif