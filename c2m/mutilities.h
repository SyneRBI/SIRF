#ifndef IUTILITIES_TO_MATLAB_INTERFACE
#define IUTILITIES_TO_MATLAB_INTERFACE

#define IUTILITIES_FOR_MATLAB
#include "shrhelp.h"

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
