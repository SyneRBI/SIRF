#define SPTR(Base, X, Object) shared_ptr< Base > X(new Object)
#define NEW_SPTR(Base, X, Object) \
	shared_ptr< Base >* X = new shared_ptr< Base >(new Object)
#define NEW_SPTR_FROM_PTR(Object, X, P) \
	shared_ptr< Object >* X = new shared_ptr< Object >(P)
#define SPTR_FROM_HANDLE(Object, X, H) \
	shared_ptr<Object> X = objectSptrFromHandle<Object>(H);

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

template<class Base, class Object>
static void*
newObjectHandle()
{
	NEW_SPTR(Base, ptr_sptr, Object);
	ObjectHandle<Base>* ptr_handle = new ObjectHandle<Base>(*ptr_sptr);
	delete ptr_sptr;
	return (void*)ptr_handle;
}

template<class Base>
static void*
newObjectHandle(shared_ptr<Base>* ptr_sptr)
{
	ObjectHandle<Base>* ptr_handle = new ObjectHandle<Base>(*ptr_sptr);
	delete ptr_sptr;
	return (void*)ptr_handle;
}

template<class T>
void*
sptrObjectHandle(shared_ptr<T> sptr) {
	ObjectHandle<T>* ptr_handle = new ObjectHandle<T>(sptr);
	return (void*)ptr_handle;
}

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

template<class Base, class Object>
Object&
objectFromHandle(const DataHandle* handle) {
	void* ptr = handle->data();
	if (ptr == 0)
		THROW("zero data pointer cannot be dereferenced");
	CAST_PTR(shared_ptr<Base>, ptr_sptr, ptr);
	if (!ptr_sptr->get())
		THROW("zero object pointer cannot be dereferenced");
	CAST_PTR(Object, ptr_object, ptr_sptr->get());
	return *ptr_object;
}

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
