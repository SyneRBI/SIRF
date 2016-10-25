#include <complex>

#include "data_handle.h"

char* charDataFromDataHandle(const DataHandle* ptr_h) 
{
	void* ptr_d = ptr_h->data();
	if (!ptr_d)
		return 0;
	else
		return (char*)ptr_d;
}

extern "C" {

	void* newDataHandle() 
	{
		return (void*)new DataHandle;
	}
	void deleteDataHandle(void* ptr) 
	{
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

	void* copyOfObjectHandle(void* ptr)
	{
		try {
			CAST_PTR(anObjectHandle, ptr_obj, ptr);
			return (void*)ptr_obj->copy();
		}
		catch (LocalisedException& le) {
			ExecutionStatus status(le.what(), le.file(), le.line());
			DataHandle* handle = new DataHandle;
			handle->set(0, &status);
			return (void*)handle;
		}
	}

	void* charDataHandle(const char* s) 
	{
		DataHandle* h = new DataHandle;
		size_t len = strlen(s);
		char* d = (char*)malloc(len + 1);
		//strcpy_s(d, len + 1, s);
		strcpy(d, s);
		h->set((void*)d, 0, GRAB);
		return (void*)h;
	}

	void* intDataHandle(int i)
	{
		return dataHandle<int>(i);
	}	
	void* floatDataHandle(float i) 
	{
		return dataHandle<float>(i);
	}
	void* doubleDataHandle(double i)
	{
		return dataHandle<double>(i);
	}

	char* charDataFromHandle(const void* ptr)
	{
		return charDataFromDataHandle((const DataHandle*)ptr);
	}
	int intDataFromHandle(const void* ptr)
	{
		return dataFromHandle<int>(ptr);
	}
	float floatDataFromHandle(const void* ptr) 
	{
		return dataFromHandle<float>(ptr);
	}
	double doubleDataFromHandle(const void* ptr)
	{
		return dataFromHandle<double>(ptr);
	}

	double doubleReDataFromHandle(const void* ptr)
	{
		std::complex<double> z = dataFromHandle<std::complex<double> >(ptr);
		return z.real();
	}
	double doubleImDataFromHandle(const void* ptr)
	{
		std::complex<double> z = dataFromHandle<std::complex<double> >(ptr);
		return z.imag();
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
