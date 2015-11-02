#include "dh.h"

extern "C" {
	int executionStatus(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		return (ptr_h->status() ? 1 : 0);
	}

	const char* executionError(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		if (ptr_h->status())
			return ptr_h->status()->error();
		else
			return "";
	}

	const char* executionErrorFile(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		if (ptr_h->status())
			return ptr_h->status()->file();
		else
			return "";
	}

	int executionErrorLine(const void* ptr) {
		const DataHandle* ptr_h = (const DataHandle*)ptr;
		if (ptr_h->status())
			return ptr_h->status()->line();
		else
			return 0;
	}
}
