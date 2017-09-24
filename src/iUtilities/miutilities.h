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
#ifndef IUTILITIES_TO_MATLAB_INTERFACE
#define IUTILITIES_TO_MATLAB_INTERFACE

#define IUTILITIES_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

#ifndef IUTILITIES_FOR_MATLAB
 extern "C" {
#endif
EXPORTED_FUNCTION  void* mNewDataHandle();
EXPORTED_FUNCTION 	void mDeleteDataHandle(void* ptr);
EXPORTED_FUNCTION 	void* mCopyOfObjectHandle(void* ptr);
EXPORTED_FUNCTION 	void mDeleteObject(void* ptr);
EXPORTED_FUNCTION 	void* mCharDataHandle(const char* s);
EXPORTED_FUNCTION 	void* mIntDataHandle(int i);
EXPORTED_FUNCTION 	void* mFloatDataHandle(float i);
EXPORTED_FUNCTION 	void* mDoubleDataHandle(double i);
EXPORTED_FUNCTION 	char* mCharDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	int mIntDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	float mFloatDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	float mFloatReDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	float mFloatImDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	double mDoubleDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	double mDoubleReDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	double mDoubleImDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	int mExecutionStatus(const void* ptr);
EXPORTED_FUNCTION 	const char* mExecutionError(const void* ptr);
EXPORTED_FUNCTION 	const char* mExecutionErrorFile(const void* ptr);
EXPORTED_FUNCTION 	int mExecutionErrorLine(const void* ptr);
#ifndef IUTILITIES_FOR_MATLAB
}
#endif

#endif
