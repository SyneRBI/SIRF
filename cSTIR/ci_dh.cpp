#include <iostream>
#include <string>
using namespace std;

#include "dh.h"

char* charDataFromHandle(const DataHandle* ptr_h) {
	void* ptr_d = ptr_h->data();
	if (!ptr_d)
		return 0;
	else
		return (char*)ptr_d;
}

extern "C" {
	void* newDataHandle() {
		return (void*)new DataHandle;
	}

	void* refDataHandle(void* ptr) {
		DataHandle* h = (DataHandle*)ptr;
		DataHandle* hr = new DataHandle;
		hr->set(h->data(), h->status());
		return (void*)hr;
	}

	void* charDataHandle(const char* s) {
		DataHandle* h = new DataHandle;
		size_t len = strlen(s);
		char* d = (char*)malloc(len + 1);
		//strcpy_s(d, len + 1, s);
		strcpy(d, s);
		h->set((void*)d, 0, GRAB);
		return (void*)h;
	}
	void* intDataHandle(int i) {
		return dataHandle<int>(i);
	}
	void* floatDataHandle(float i) {
		return dataHandle<float>(i);
	}
	void* doubleDataHandle(double i) {
		return dataHandle<double>(i);
	}

	char* charDataFromHandle(const void* ptr) {
		return charDataFromHandle((const DataHandle*)ptr);
	}
	int intDataFromHandle(const void* ptr) {
		return dataFromHandle<int>(ptr);
	}
	float floatDataFromHandle(const void* ptr) {
		return dataFromHandle<float>(ptr);
	}
	double doubleDataFromHandle(const void* ptr) {
		return dataFromHandle<double>(ptr);
	}

	void deleteDataHandle(void* ptr) {
		if (ptr)
			delete (DataHandle*)ptr;
	}
}
