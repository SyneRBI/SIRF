/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

/*!
\file
\ingroup C Interface to C++ Objects
\brief Execution status type and basic wrapper for C++ objects.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef DATA_HANDLE_TYPES
#define DATA_HANDLE_TYPES

#include <stdlib.h>
#include <string>

#include "localised_exception.h"

#define NEW(T, X) T* X = new T
#define CAST_PTR(T, X, Y) T* X = (T*)Y
#define THROW(msg) throw LocalisedException(msg, __FILE__, __LINE__)
#define CATCH \
	catch (LocalisedException& se) {\
		ExecutionStatus status(se);\
		DataHandle* handle = new DataHandle;\
		handle->set(0, &status);\
		return (void*)handle;\
				}\
	catch (std::string msg) {\
		ExecutionStatus status(msg.c_str(), __FILE__, __LINE__);\
		DataHandle* handle = new DataHandle;\
		handle->set(0, &status);\
		return (void*)handle;\
        }\
	catch (...) {\
		ExecutionStatus status("unhandled", __FILE__, __LINE__);\
		DataHandle* handle = new DataHandle;\
		handle->set(0, &status);\
		return (void*)handle;\
				}\

/*!
\ingroup C Interface to C++ Objects
\brief Execution status type.

An ExecutionStatus object is created when an exception is caught (see above).
It stores the exeption's error message and position (file name and line number). 
*/
class ExecutionStatus {
public:
	ExecutionStatus() : _error(0), _file(0), _line(0) {}
	ExecutionStatus(const char* error, const char* file, int line) {
		set(error, file, line);
	}
	ExecutionStatus(const ExecutionStatus& s) {
		set(s.error(), s.file(), s.line());
	}
	ExecutionStatus(const LocalisedException& ex) {
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

/*!
\ingroup C Interface to C++ Objects
\brief Basic wrapper for C++ objects.

A DataHandle object stores data address (void* _data) and the current
execution status (ExecutionStatus _status).
SIRF C interface functions work with pointers to DataHandle objects
cast to void*.
*/
class DataHandle {
public:
	DataHandle() : _data(0), _status(0), _owns_data(0) {}
	virtual ~DataHandle() {
		if (_data && _owns_data)
			free(_data);
		delete _status;
	}
	void set(void* data, const ExecutionStatus* status = 0, int grab = 0) {
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
	bool _owns_data; // can free _data
	void* _data; // data address
	ExecutionStatus* _status; // execution status
};

#define GRAB 1

/*!
\ingroup C Interface to C++ Objects
\brief Data wrapper.

Wraps an object of type T into DataHandle.
The data is owned by the DataHandle object and hence will be deleted by its 
destructor.
*/
template <typename T>
void
setDataHandle(DataHandle* h, T x)
{
	T* ptr = (T*)malloc(sizeof(T));
	*ptr = x;
	h->set((void*)ptr, 0, GRAB);
}

/*!
\ingroup C Interface to C++ Objects
\brief Data wrapper constructor.

Creates a new DataHandle to wrap an object of type T. 
*/
template <typename T>
void*
dataHandle(T x)
{
	DataHandle* h = new DataHandle;
	setDataHandle<T>(h, x);
	return (void*)h;
}

/*!
\ingroup C Interface to C++ Objects
\brief Data extractor.

Returns a copy of the data stored in a DataHandle object.
*/
template <typename T>
T // must have a proper copy constructor
dataFromHandle(const void* ptr)
{
	DataHandle* ptr_h = (DataHandle*)ptr;
	void* ptr_d = ptr_h->data();
	if (!ptr_d)
		return 0;
	else
		return *((T*)ptr_d);
}

// yet another kludge to stop matlab on linux from crashing

inline char* charDataFromDataHandle(const DataHandle* ptr_h)
{
	void* ptr_d = ptr_h->data();
	if (!ptr_d)
		return 0;
	else
		return (char*)ptr_d;
}

inline void* charDataHandleFromCharData(const char* s)
{
	DataHandle* h = new DataHandle;
	size_t len = strlen(s);
	char* d = (char*)malloc(len + 1);
	//strcpy_s(d, len + 1, s);
	strcpy(d, s);
	h->set((void*)d, 0, GRAB);
	return (void*)h;
}

//#define SPTR_NAMESPACE std_sptr
//#include "object_handle.h"

#endif
