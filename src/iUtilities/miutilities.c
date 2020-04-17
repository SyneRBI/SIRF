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
#define IUTILITIES_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

#include <mex.h>
#include "matrix.h"
#include "sirf/iUtilities/iutilities.h"

#ifndef IUTILITIES_FOR_MATLAB
 extern "C" {
#endif
EXPORTED_FUNCTION  void* mNewDataHandle() {
	return newDataHandle();
}
EXPORTED_FUNCTION 	void mDeleteDataHandle(void* ptr) {
	deleteDataHandle(ptr);
}
EXPORTED_FUNCTION 	void* mCharDataHandle(const char* s) {
	return charDataHandle(s);
}
EXPORTED_FUNCTION 	void* mIntDataHandle(int i) {
	return intDataHandle(i);
}
EXPORTED_FUNCTION 	void* mFloatDataHandle(float i) {
	return floatDataHandle(i);
}
EXPORTED_FUNCTION 	void* mDoubleDataHandle(double i) {
	return doubleDataHandle(i);
}
EXPORTED_FUNCTION 	char* mCharDataFromHandle(const void* ptr) {
	return charDataFromHandle(ptr);
}
EXPORTED_FUNCTION     bool mBoolDataFromHandle(const void* ptr) {
	return boolDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	int mIntDataFromHandle(const void* ptr) {
	return intDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	int mIntDataItemFromHandle(const void* ptr, int i) {
	return intDataItemFromHandle(ptr, i);
}
EXPORTED_FUNCTION 	int mUint16DataItemFromHandle(const void* ptr, int i) {
	return uint16DataItemFromHandle(ptr, i);
}
EXPORTED_FUNCTION 	int mUint32DataItemFromHandle(const void* ptr, int i) {
	return uint32DataItemFromHandle(ptr, i);
}
EXPORTED_FUNCTION 	int mUint64DataItemFromHandle(const void* ptr, int i) {
	return uint64DataItemFromHandle(ptr, i);
}
EXPORTED_FUNCTION 	float mFloatDataFromHandle(const void* ptr) {
	return floatDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	float mFloatDataItemFromHandle(const void* ptr, int i) {
	return floatDataItemFromHandle(ptr, i);
}
EXPORTED_FUNCTION 	float mFloatReDataFromHandle(const void* ptr) {
	return floatReDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	float mFloatImDataFromHandle(const void* ptr) {
	return floatImDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	double mDoubleDataFromHandle(const void* ptr) {
	return doubleDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	double mDoubleReDataFromHandle(const void* ptr) {
	return doubleReDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	double mDoubleImDataFromHandle(const void* ptr) {
	return doubleImDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	int mExecutionStatus(const void* ptr) {
	return executionStatus(ptr);
}
EXPORTED_FUNCTION 	const char* mExecutionError(const void* ptr) {
	return executionError(ptr);
}
EXPORTED_FUNCTION 	const char* mExecutionErrorFile(const void* ptr) {
	return executionErrorFile(ptr);
}
EXPORTED_FUNCTION 	int mExecutionErrorLine(const void* ptr) {
	return executionErrorLine(ptr);
}
#ifndef IUTILITIES_FOR_MATLAB
}
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {}
