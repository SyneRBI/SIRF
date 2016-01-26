#ifndef CGADGETRON_TO_MATLAB_INTERFACE
#define CGADGETRON_TO_MATLAB_INTERFACE

#define CGADGETRON_FOR_MATLAB
#include "shrhelp.h"

EXPORTED_FUNCTION void* mCopyOfObject(void* ptr);
EXPORTED_FUNCTION void mDeleteObject(void* ptr);
EXPORTED_FUNCTION void* mNewDataHandle();
EXPORTED_FUNCTION void mDeleteDataHandle(void* ptr);
EXPORTED_FUNCTION int mExecutionStatus(const void* ptr);
EXPORTED_FUNCTION const char* mExecutionError(const void* ptr);
EXPORTED_FUNCTION const char* mExecutionErrorFile(const void* ptr);
EXPORTED_FUNCTION int mExecutionErrorLine(const void* ptr);
EXPORTED_FUNCTION void* mNewTextPrinter(const char* stream);
EXPORTED_FUNCTION void* mNewTextWriter(const char* stream);
EXPORTED_FUNCTION void mOpenChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION void mCloseChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION void mSetWriter(void* ptr_w, int channel);
EXPORTED_FUNCTION void mResetWriter();
EXPORTED_FUNCTION void mPrintText(const char* text);
EXPORTED_FUNCTION void mDeleteTextPrinter(void* ptr);
EXPORTED_FUNCTION void mDeleteTextWriter(void* ptr_w);

EXPORTED_FUNCTION void* mNewMexPrinter();
EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr);
#endif
