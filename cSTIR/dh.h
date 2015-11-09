#ifndef DATA_HANDLE
#define DATA_HANDLE

#include <stdlib.h>

#include <boost\algorithm\string.hpp>

#include "StirException.h"

class ExecutionStatus {
public:
	ExecutionStatus() : _error(0), _file(0), _line(0) {}
	ExecutionStatus(const char* error, const char* file, int line) {
		set(error, file, line);
	}
	ExecutionStatus(const ExecutionStatus& s) {
		set(s.error(), s.file(), s.line());
	}
	ExecutionStatus(const StirException& ex) {
		set(ex.what(), ex.file(), ex.line());
	}
	~ExecutionStatus() {
		delete[] _error;
		delete[] _file;
	}
	const char* error() const { return _error; }
	const char* file() const { return _file; }
	int line() const { return _line; }
private:
	char* _error;
	char* _file;
	int _line;
	void set(const char* error, const char* file, int line) {
		size_t size;
		if (error) {
			size = strlen(error) + 1;
			_error = new char[size];
			memcpy(_error, error, size);
		}
		else
			_error = 0;
		if (file) {
			size = strlen(file) + 1;
			_file = new char[size];
			memcpy(_file, file, size);
		}
		else
			_file = 0;
		_line = line;
	}
};

class DataHandle {
public:
	DataHandle() : _data(0), _status(0), _owns_data(0) {}
	virtual ~DataHandle() {
		if (_data && _owns_data)
			free(_data);
		delete _status;
	}
	void set(void* data, const ExecutionStatus* status, int grab = 0) {
		if (status) {
			delete _status;
			_status = new ExecutionStatus(*status);
		}
		if (_data && _owns_data)
			free(_data);
		_data = data;
		_owns_data = grab != 0;
	}
	void* data() const { return _data; }
	const ExecutionStatus* status() const { return _status; }
protected:
	bool _owns_data;
	void* _data;
	ExecutionStatus* _status;
};

#define GRAB 1

template <typename T>
void
setDataHandle(DataHandle* h, T x)
{
	T* ptr = (T*)malloc(sizeof(T));
	*ptr = x;
	h->set((void*)ptr, 0, GRAB);
}

template <typename T>
void*
dataHandle(T x)
{
	DataHandle* h = new DataHandle;
	setDataHandle<T>(h, x);
	return (void*)h;
}

template <typename T>
T
dataFromHandle(const void* ptr)
{
	DataHandle* ptr_h = (DataHandle*)ptr;
	void* ptr_d = ptr_h->data();
	if (!ptr_d)
		return 0;
	else
		return *((T*)ptr_d);
}

#endif