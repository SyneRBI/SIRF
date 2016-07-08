#include <mex.h>
#define EXPORT_FCNS
#define IUTILITIES_FOR_MATLAB
#include "matrix.h"
#include "shrhelp.h"
#include "iutilities.h"

#ifndef IUTILITIES_FOR_MATLAB
extern "C" {
#endif
EXPORTED_FUNCTION  void* mNewDataHandle() {
	return newDataHandle();
}
EXPORTED_FUNCTION 	void mDeleteDataHandle(void* ptr) {
	deleteDataHandle(ptr);
}
EXPORTED_FUNCTION 	void* mCopyOfObjectHandle(void* ptr) {
	return copyOfObjectHandle(ptr);
}
EXPORTED_FUNCTION 	void mDeleteObject(void* ptr) {
	deleteObject(ptr);
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
EXPORTED_FUNCTION 	int mIntDataFromHandle(const void* ptr) {
	return intDataFromHandle(ptr);
}
EXPORTED_FUNCTION 	float mFloatDataFromHandle(const void* ptr) {
	return floatDataFromHandle(ptr);
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
