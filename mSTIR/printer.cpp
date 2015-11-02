#include "printer.h"

extern "C" {
	void* newMexPrinter() {
		return (void*)new mexTextPrinter;
	}
	void deleteMexPrinter(void* ptr) {
		delete (mexTextPrinter*)ptr;
	}
}

