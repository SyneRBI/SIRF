#ifndef MEX_PRINTER
#define MEX_PRINTER

#include <mex.h>

#include "stir/TextWriter.h"

USING_NAMESPACE_STIR

class mexTextPrinter : public aTextWriter {
public:
	virtual void write(const char* text) const {
        //mexPrintf("mexPrintf is called...\n");
		mexPrintf(text);
	}
};

#endif
