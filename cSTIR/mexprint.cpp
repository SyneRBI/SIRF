#include "mexprint.h"

extern "C" {
	void* newMexPrinter() {
		//        mexPrintf("creating mexTextPrinter\n");
		return (void*)new mexTextPrinter;
	}
	void deleteMexPrinter(void* ptr) {
		//        mexPrintf("deleting mexTextPrinter...");
		delete (mexTextPrinter*)ptr;
		//        mexPrintf("ok\n");
	}
}

