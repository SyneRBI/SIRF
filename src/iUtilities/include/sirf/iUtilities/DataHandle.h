/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2018 - 2019 University College London
This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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
\author SyneRBI
*/

#ifndef DATA_HANDLE_TYPES
#define DATA_HANDLE_TYPES

#include <stdlib.h>
#include <string>
#include <vector>

#include "sirf/iUtilities/LocalisedException.h"

#define NEW(T, X) T* X = new T
#define CAST_PTR(T, X, Y) T* X = (T*)Y
//#define THROW(msg) throw LocalisedException(msg, __FILE__, __LINE__)
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
    catch (const std::exception &error) {\
        ExecutionStatus status(error.what(), __FILE__, __LINE__);\
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

/// Typedef of vector of void pointers for a vector of handles
typedef std::vector<void const *> DataHandleVector;

/*!
\ingroup C Interface to C++ Objects
\brief Execution status type.

An ExecutionStatus object is created when an exception is caught (see above).
It stores the exeption's error message and position (file name and line number). 
*/
class ExecutionStatus {
public:
	ExecutionStatus() : _line(0) {}
	ExecutionStatus(const std::string& error, const std::string& file, int line) {
		set(error, file, line);
	}
	ExecutionStatus(const ExecutionStatus& s) {
		set(s.error(), s.file(), s.line());
	}
	ExecutionStatus(const LocalisedException& ex) {
		set(ex.what(), ex.file(), ex.line());
	}
    ~ExecutionStatus() {}
	const std::string& error() const { return _error; }
	const std::string& file() const { return _file; }
	int line() const { return _line; }
private:
	std::string _error;
	std::string _file;
	int _line;
	void set(const std::string& error, const std::string& file, int line) {
        _error = error;
        _file = file;
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
	void set_status(const std::string& error, const std::string& file, int line)
	{
		if (_status)
			delete _status;
		_status = new ExecutionStatus(error, file, line);
	}
	void* data() const { return _data; }
	const ExecutionStatus* status() const { return _status; }
protected:
	bool _owns_data; // can free _data
	void* _data; // data address
	ExecutionStatus* _status; // execution status
};

#include <boost/shared_ptr.hpp>

template<class Base>
class ObjectHandle : public DataHandle {
public:
	ObjectHandle(const ObjectHandle& obj) {
		if (obj.uses_boost_sptr()) {
			NEW(boost::shared_ptr<Base>, ptr_sptr);
			*ptr_sptr = *(boost::shared_ptr<Base>*)obj.data();
			_data = (void*)ptr_sptr;
		}
		else {
			NEW(std::shared_ptr<Base>, ptr_sptr);
			*ptr_sptr = *(std::shared_ptr<Base>*)obj.data();
			_data = (void*)ptr_sptr;
		}
		if (obj._status)
			_status = new ExecutionStatus(*obj._status);
		else
			_status = 0;
	}
	ObjectHandle(const std::shared_ptr<Base>& sptr,
		const ExecutionStatus* status = 0) : _boost_sptr(false) {
		NEW(std::shared_ptr<Base>, ptr_sptr);
		*ptr_sptr = sptr;
		_data = (void*)ptr_sptr;
		if (status)
			_status = new ExecutionStatus(*status);
		else
			_status = 0;
	}
	ObjectHandle(const boost::shared_ptr<Base>& sptr,
		const ExecutionStatus* status = 0) : _boost_sptr(true) {
		NEW(boost::shared_ptr<Base>, ptr_sptr);
		*ptr_sptr = sptr;
		_data = (void*)ptr_sptr;
		if (status)
			_status = new ExecutionStatus(*status);
		else
			_status = 0;
	}
	virtual ~ObjectHandle() {
		delete _status;
		_status = 0;
		if (_boost_sptr) {
			CAST_PTR(boost::shared_ptr<Base>, ptr_sptr, _data);
			delete ptr_sptr;
		}
		else {
			CAST_PTR(std::shared_ptr<Base>, ptr_sptr, _data);
			delete ptr_sptr;
		}
	}
	bool uses_boost_sptr() const
	{
		return _boost_sptr;
	}
protected:
	bool _boost_sptr;
};

template<class Object>
static void*
newObjectHandle(std::shared_ptr<Object> sptr)
{
	return (void*)new ObjectHandle<Object>(sptr);
}

template<class Object>
static void*
newObjectHandle(boost::shared_ptr<Object> sptr)
{
	return (void*)new ObjectHandle<Object>(sptr);
}

template<class Object>
Object&
objectFromHandle(const void* h) {
	ObjectHandle<Object>* handle = (ObjectHandle<Object>*)h;
	void* ptr = handle->data();
	if (ptr == 0)
		THROW("zero data pointer cannot be dereferenced");
	if (handle->uses_boost_sptr()) {
		CAST_PTR(boost::shared_ptr<Object>, ptr_sptr, ptr);
		Object* ptr_obj = ptr_sptr->get();
		if (ptr_obj == 0)
			THROW("zero object pointer cannot be dereferenced");
		return *ptr_obj;
	}
	else {
		CAST_PTR(std::shared_ptr<Object>, ptr_sptr, ptr);
		Object* ptr_obj = ptr_sptr->get();
		if (ptr_obj == 0)
			THROW("zero object pointer cannot be dereferenced");
		return *ptr_obj;
	}
}

template<class Object>
void
getObjectSptrFromHandle(const void* h, std::shared_ptr<Object>& sptr) {
	ObjectHandle<Object>* handle = (ObjectHandle<Object>*)h;
	if (handle->uses_boost_sptr())
		THROW("cannot cast boost::shared_ptr to std::shared_ptr");
	void* ptr = handle->data();
	if (ptr == 0)
		THROW("zero data pointer cannot be dereferenced");
	CAST_PTR(std::shared_ptr<Object>, ptr_sptr, ptr);
	sptr = *ptr_sptr;
}

template<class Object>
void
getObjectSptrFromHandle(const void* h, boost::shared_ptr<Object>& sptr) {
	ObjectHandle<Object>* handle = (ObjectHandle<Object>*)h;
	if (!handle->uses_boost_sptr())
		THROW("cannot cast std::shared_ptr to boost::shared_ptr");
	void* ptr = handle->data();
	if (ptr == 0)
		THROW("zero data pointer cannot be dereferenced");
	CAST_PTR(boost::shared_ptr<Object>, ptr_sptr, ptr);
	sptr = *ptr_sptr;
}

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

#endif
