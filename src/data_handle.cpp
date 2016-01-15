#include "data_handle.h"

extern "C" {

	void* newDataHandle() {
		return (void*)new DataHandle;
	}
	
	void deleteDataHandle(void* ptr) {
		if (ptr)
			delete (DataHandle*)ptr;
	}
	
	void deleteObject(void* ptr)
	{
		if (!ptr)
			return;
		CAST_PTR(anObjectHandle, ptr_obj, ptr);
		delete ptr_obj;
	}

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