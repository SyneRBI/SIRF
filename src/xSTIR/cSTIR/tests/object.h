#pragma once

//#define CAST_REF(T, X, Y) T& X = (T&)Y
#define CREATE_OBJ(Obj, X, sptr_X, Par) \
	boost::shared_ptr< Obj > sptr_X(new Obj(Par)); \
	Obj& X = (Obj&)*sptr_X
#define CREATE_OBJECT(Base, Object, X, sptr_X, Par) \
	boost::shared_ptr< Base > sptr_X(new Object(Par)); \
	Object& X = (Object&)*sptr_X
