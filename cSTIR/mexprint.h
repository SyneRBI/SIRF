#ifndef MEX_PRINTER
#define MEX_PRINTER

#include <mex.h>

#include "TextWriter.h"

class mexTextPrinter : public aTextWriter {
public:
	virtual void write(const char* text) const {
		//mexPrintf("mexPrintf is called...\n");
		mexPrintf(text);
	}
};

#endif
