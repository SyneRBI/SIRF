#include <mex.h>
#define EXPORT_FCNS
#define CGADGETRON_FOR_MATLAB
#include "matrix.h"
#include "shrhelp.h"
#include "iutilities.h"

#ifndef CGADGETRON_FOR_MATLAB
extern "C" {
#endif
EXPORTED_FUNCTION  void* mNewDataHandle() {
	return newDataHandle();
}
EXPORTED_FUNCTION 	void* mCopyOfObjectHandle(void* ptr) {
	return copyOfObjectHandle(ptr);
}
EXPORTED_FUNCTION 	void mDeleteObject(void* ptr) {
	deleteObject(ptr);
}
EXPORTED_FUNCTION 	void mDeleteDataHandle(void* ptr) {
	deleteDataHandle(ptr);
}
EXPORTED_FUNCTION 	void* mIntDataHandle(int i) {
	return intDataHandle(i);
}
EXPORTED_FUNCTION 	void* mDoubleDataHandle(double i) {
	return doubleDataHandle(i);
}
EXPORTED_FUNCTION 	int mIntDataFromHandle(const void* ptr) {
	return intDataFromHandle(ptr);
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
EXPORTED_FUNCTION 	void* mNewTextPrinter(const char* stream) {
	return newTextPrinter(stream);
}
EXPORTED_FUNCTION 	void* mNewTextWriter(const char* stream) {
	return newTextWriter(stream);
}
EXPORTED_FUNCTION 	void mOpenChannel(int channel, void* ptr_w) {
	openChannel(channel, ptr_w);
}
EXPORTED_FUNCTION 	void mCloseChannel(int channel, void* ptr_w) {
	closeChannel(channel, ptr_w);
}
EXPORTED_FUNCTION 	void mSetWriter(void* ptr_w, int channel) {
	setWriter(ptr_w, channel);
}
EXPORTED_FUNCTION 	void mResetWriter() {
	resetWriter();
}
EXPORTED_FUNCTION 	void mPrintText(const char* text) {
	printText(text);
}
EXPORTED_FUNCTION 	void mDeleteTextPrinter(void* ptr) {
	deleteTextPrinter(ptr);
}
EXPORTED_FUNCTION 	void mDeleteTextWriter(void* ptr_w) {
	deleteTextWriter(ptr_w);
}
#ifndef CGADGETRON_FOR_MATLAB
}
#endif


void* newMexPrinter();
void deleteMexPrinter(void* ptr);
EXPORTED_FUNCTION void* mNewMexPrinter() {
  return newMexPrinter();
}
EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr) {
  deleteMexPrinter(ptr);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {}
