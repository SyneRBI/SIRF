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
\brief A wrapper to a shared pointer.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#define SPTR(Base, X, Object) shared_ptr< Base > X(new Object)
#define NEW_SPTR(Base, X, Object) \
	shared_ptr< Base >* X = new shared_ptr< Base >(new Object)
#define NEW_SPTR_FROM_PTR(Object, X, P) \
	shared_ptr< Object >* X = new shared_ptr< Object >(P)
#define SPTR_FROM_HANDLE(Object, X, H) \
	shared_ptr<Object> X = objectSptrFromHandle<Object>(H);

/*!
\ingroup C Interface to C++ Objects
\brief A wrapper to a shared pointer.

A derived DataHandle type that stores a shared pointer, _data being 
a pointer to it.
*/
template<class Base>
class ObjectHandle : public DataHandle {
public:
	ObjectHandle(const ObjectHandle& obj) {
		NEW(shared_ptr<Base>, ptr_sptr);
		*ptr_sptr = *(shared_ptr<Base>*)obj.data();
		_data = (void*)ptr_sptr;
		if (obj._status)
			_status = new ExecutionStatus(*obj._status);
		else
			_status = 0;
	}
	ObjectHandle(const shared_ptr<Base>& sptr,
		const ExecutionStatus* status = 0) {
		NEW(shared_ptr<Base>, ptr_sptr);
		*ptr_sptr = sptr;
		_data = (void*)ptr_sptr;
		if (status)
			_status = new ExecutionStatus(*status);
		else
			_status = 0;
	}
	virtual ~ObjectHandle() {
		CAST_PTR(shared_ptr<Base>, ptr_sptr, _data);
		delete _status;
		_status = 0;
		delete ptr_sptr;
	}
};

/*!
\ingroup C Interface to C++ Objects
\brief A new wrapper for a new object.

Creates an object of the class Object derived from the class Base and 
an ObjectHandle storing a shared pointer to it.
Returns a void pointer to the ObjectHandle object.
*/
//template<class Base, class Object>
template<class Object>
static void*
newObjectHandle()
{
	//NEW_SPTR(Base, ptr_sptr, Object);
	//ObjectHandle<Base>* ptr_handle = new ObjectHandle<Base>(*ptr_sptr);
	//delete ptr_sptr;
	//return (void*)ptr_handle;
	shared_ptr<Object> sptr(new Object);
	return new ObjectHandle<Object>(sptr);
}

/*!
\ingroup C Interface to C++ Objects
\brief A new wrapper for an existing shared pointer.

Creates an ObjectHandle that stores a shared pointer pointed to by the argument.
Returns a void pointer to the ObjectHandle object.
*/
template<class Object>
static void*
newObjectHandle(shared_ptr<Object>* ptr_sptr)
{
	ObjectHandle<Object>* ptr_handle = new ObjectHandle<Object>(*ptr_sptr);
	delete ptr_sptr;
	return (void*)ptr_handle;
}

/*!
\ingroup C Interface to C++ Objects
\brief A new wrapper for an existing shared pointer.

Creates an ObjectHandle that stores a shared pointer supplied as the argument.
Returns a void pointer to the ObjectHandle object.
*/
template<class T>
void*
sptrObjectHandle(shared_ptr<T> sptr) {
	ObjectHandle<T>* ptr_handle = new ObjectHandle<T>(sptr);
	return (void*)ptr_handle;
}

/*!
\ingroup C Interface to C++ Objects
\brief An object accessor.

Returns a reference to the object wrapped up by the ObjectHandle pointed to
by the argument.
*/
template<class Object>
Object&
objectFromHandle(const void* h) {
	DataHandle* handle = (DataHandle*)h;
	void* ptr = handle->data();
	if (ptr == 0)
		THROW("zero data pointer cannot be dereferenced");
	CAST_PTR(shared_ptr<Object>, ptr_sptr, ptr);
	if (!ptr_sptr->get())
		THROW("zero object pointer cannot be dereferenced");
	CAST_PTR(Object, ptr_object, ptr_sptr->get());
	return *ptr_object;
}

/*!
\ingroup C Interface to C++ Objects
\brief An object accessor.

Returns a reference to the shared pointer to the object wrapped up by 
the ObjectHandle pointed to by the argument.
*/
template<class Object>
shared_ptr<Object>&
objectSptrFromHandle(const void* h) {
	DataHandle* handle = (DataHandle*)h;
	void* ptr = handle->data();
	if (ptr == 0)
		THROW("zero data pointer cannot be dereferenced");
	CAST_PTR(shared_ptr<Object>, ptr_sptr, ptr);
	if (!ptr_sptr->get())
		THROW("zero object pointer cannot be dereferenced");
	return *ptr_sptr;
}

template<class Base>
Base&
objectFromHandle(const DataHandle* handle) {
	void* ptr = handle->data();
	if (ptr == 0)
		THROW("zero data pointer cannot be dereferenced");
	CAST_PTR(shared_ptr<Base>, ptr_sptr, ptr);
	if (!ptr_sptr->get())
		THROW("zero object pointer cannot be dereferenced");
	CAST_PTR(Base, ptr_object, ptr_sptr->get());
	return *ptr_object;
}

//template<class Base, class Object>
//Object&
//objectFromHandle(const DataHandle* handle) {
//	void* ptr = handle->data();
//	if (ptr == 0)
//		THROW("zero data pointer cannot be dereferenced");
//	CAST_PTR(shared_ptr<Base>, ptr_sptr, ptr);
//	if (!ptr_sptr->get())
//		THROW("zero object pointer cannot be dereferenced");
//	CAST_PTR(Object, ptr_object, ptr_sptr->get());
//	return *ptr_object;
//}

template<class Base>
shared_ptr<Base>&
objectSptrFromHandle(const DataHandle* handle) {
	void* ptr = handle->data();
	if (ptr == 0)
		THROW("zero data pointer cannot be dereferenced");
	CAST_PTR(shared_ptr<Base>, ptr_sptr, ptr);
	if (!ptr_sptr->get())
		THROW("zero object pointer cannot be dereferenced");
	return *ptr_sptr;
}

template<class T>
shared_ptr<T>
sptrDataFromHandle(const DataHandle* handle) {
	return *(shared_ptr<T>*)handle->data();
}
