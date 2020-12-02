/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 Rutherford Appleton Laboratory STFC
Copyright 2020 University College London

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

#pragma once

// in case #pragma once not supported
#ifndef SIRF_NUMBER_REFERENCE
#define SIRF_NUMBER_REFERENCE

#include <complex>
#include <typeinfo>
#include <stdexcept>
/// Number type. Taken from ismrmrd/ismrmrd.h (saves having to include 
/// it in this file, which should be independent of it)
/// class NumberType
class NumberType
{
public:
	enum Type {
		USHORT = 1, /**< corresponds to uint16_t */
		SHORT = 2, /**< corresponds to int16_t */
		UINT = 3, /**< corresponds to uint32_t */
		INT = 4, /**< corresponds to int32_t */
		FLOAT = 5, /**< corresponds to float */
		DOUBLE = 6, /**< corresponds to double */
		CXFLOAT = 7, /**< corresponds to complex float */
		CXDOUBLE = 8  /**< corresponds to complex double */
	};
};
typedef std::complex<float>  complex_float_t;
typedef std::complex<double> complex_double_t;

namespace sirf {

	class ANumRef {
	public:
		virtual complex_double_t complex_double() const = 0;
		virtual complex_float_t complex_float() const = 0;
		virtual operator float() const = 0;
		virtual void assign(const ANumRef& ref) = 0;
		ANumRef& operator=(const ANumRef& ref)
		{
			assign(ref);
			return *this;
		}
		virtual void set_ptr(void* ptr) = 0;
        virtual NumberType::Type get_typeID() const = 0;
	};

	class FloatRef : public ANumRef {
	public:
		FloatRef(float* ptr = 0, int dummy = 0) : ptr_(ptr)
		{}
		virtual complex_double_t complex_double() const
		{
			return complex_double_t(*ptr_);
		}
		virtual complex_float_t complex_float() const
		{
			return complex_float_t(*ptr_);
		}
		virtual operator float() const
		{
			return *ptr_;
		}
		template<typename T>
		FloatRef& operator=(T v)
		{
			*ptr_ = v;
			return *this;
		}
		virtual void assign(const ANumRef& a_ref)
		{
			const FloatRef& ref = (const FloatRef&)a_ref;
			*ptr_ = float(ref);
		}
		ANumRef& operator=(const ANumRef& ref)
		{
			assign(ref);
			return *this;
		}
		void set_ptr(void* ptr)
		{
			ptr_ = (float*)ptr;
		}
		void copy(const FloatRef& ref)
		{
			ptr_ = ref.ptr_;
		}
        virtual NumberType::Type get_typeID() const
        {
            return NumberType::FLOAT;
        }
	private:
		float* ptr_;
	};

	template <typename Type>
	NumberType::Type TypeID(Type t)
	{
		if (typeid(Type) == typeid(complex_double_t))
			return NumberType::CXDOUBLE;
		else if (typeid(Type) == typeid(complex_float_t))
			return NumberType::CXFLOAT;
		else if (typeid(Type) == typeid(double))
			return NumberType::DOUBLE;
		else if (typeid(Type) == typeid(float))
			return NumberType::FLOAT;
		else if (typeid(Type) == typeid(int))
			return NumberType::INT;
		else if (typeid(Type) == typeid(unsigned int))
			return NumberType::UINT;
		else if (typeid(Type) == typeid(short))
			return NumberType::SHORT;
		else if (typeid(Type) == typeid(unsigned short))
			return NumberType::USHORT;
		else
			throw std::invalid_argument
			(std::string("unsupported numeric type ") + typeid(Type).name());
	}

	class NumRef : public ANumRef {
	public:
		NumRef(void* ptr = 0, int type = (int)NumberType::FLOAT) :
			ptr_(ptr), abs_(true), type_(type)
		{}
		NumRef(const NumRef& ref) :
			ptr_(ref.ptr_), abs_(ref.abs_), type_(ref.type_)
		{}
		void set_complex_to_real_mode(char m)
		{
			abs_ = (m == 'a');
		}
		virtual void set_ptr(void* ptr)
		{
			ptr_ = ptr;
		}
		virtual void copy(const NumRef& ref)
		{
			ptr_ = ref.ptr_;
			abs_ = ref.abs_;
			type_ = ref.type_;
		}
		complex_double_t complex_double() const
		{
			complex_double_t z;
			switch (type_) {
			case NumberType::CXDOUBLE:
				z = *(complex_double_t*)ptr_;
				break;
			case NumberType::CXFLOAT:
				z = *(complex_float_t*)ptr_;
				break;
			case NumberType::DOUBLE:
				z = float(*(double*)ptr_);
				break;
			case NumberType::FLOAT:
				z = *(float*)ptr_;
				break;
			case NumberType::INT:
				z = complex_float_t(float(*(int*)ptr_));
				break;
			case NumberType::UINT:
				z = complex_float_t(float(*(unsigned int*)ptr_));
				break;
			case NumberType::SHORT:
				z = complex_float_t(float(*(short*)ptr_));
				break;
			case NumberType::USHORT:
				z = complex_float_t(float(*(unsigned short*)ptr_));
			}
			return z;
		}
		virtual complex_float_t complex_float() const
		{
			complex_float_t z;
			switch (type_) {
			case NumberType::CXDOUBLE:
				z = *(complex_double_t*)ptr_;
				break;
			case NumberType::CXFLOAT:
				z = *(complex_float_t*)ptr_;
				break;
			case NumberType::DOUBLE:
				z = float(*(double*)ptr_);
				break;
			case NumberType::FLOAT:
				z = *(float*)ptr_;
				break;
			case NumberType::INT:
				z = complex_float_t(float(*(int*)ptr_));
				break;
			case NumberType::UINT:
				z = complex_float_t(float(*(unsigned int*)ptr_));
				break;
			case NumberType::SHORT:
				z = complex_float_t(float(*(short*)ptr_));
				break;
			case NumberType::USHORT:
				z = complex_float_t(float(*(unsigned short*)ptr_));
			}
			return z;
		}
		virtual operator float() const
		{
			float v;
			complex_float_t c;
			complex_double_t z;
			switch (type_) {
			case NumberType::CXDOUBLE:
				z = *(complex_double_t*)ptr_;
				v = float(abs_ ? abs(z) : z.real());
				break;
			case NumberType::CXFLOAT:
				c = *(complex_float_t*)ptr_;
				v = abs_ ? abs(c) : c.real();
				break;
			case NumberType::DOUBLE:
				v = float(*(double*)ptr_);
				break;
			case NumberType::FLOAT:
				v = *(float*)ptr_;
				break;
			case NumberType::INT:
				v = float(*(int*)ptr_);
				break;
			case NumberType::UINT:
				v = float(*(unsigned int*)ptr_);
				break;
			case NumberType::SHORT:
				v = float(*(short*)ptr_);
				break;
			case NumberType::USHORT:
				v = float(*(unsigned short*)ptr_);
			}
			return v;
		}
		NumRef& operator=(const NumRef& ref)
		{
			assign(ref);
			return *this;
		}
		virtual void assign(const ANumRef& a_ref)
		{
			const NumRef& ref = (const NumRef&)a_ref;
			switch (type_) {
			case NumberType::CXDOUBLE:
				*(complex_double_t*)ptr_ = ref.complex_double();
				break;
			case NumberType::CXFLOAT:
				*(complex_float_t*)ptr_ = ref.complex_float();
				break;
			case NumberType::DOUBLE:
				*(double*)ptr_ = double(ref);
				break;
			case NumberType::FLOAT:
				*(float*)ptr_ = float(ref);
				break;
			case NumberType::INT:
				*(int*)ptr_ = int(ref);
				break;
			case NumberType::UINT:
				*(unsigned int*)ptr_ = (unsigned int)ref;
				break;
			case NumberType::SHORT:
				*(short*)ptr_ = short(ref);
				break;
			case NumberType::USHORT:
				*(unsigned short*)ptr_ = (unsigned short)ref;
			}
		}
		template <typename T>
		NumRef& operator=(std::complex<T> v)
		{
			switch (type_) {
			case NumberType::CXDOUBLE:
				*(complex_double_t*)ptr_ = complex_double_t(v);
				break;
			case NumberType::CXFLOAT:
				*(complex_float_t*)ptr_ = complex_float_t(v);
				break;
			case NumberType::DOUBLE:
				*(double*)ptr_ = double(abs_ ? abs(v) : v.real());
				break;
			case NumberType::FLOAT:
				*(float*)ptr_ = float(abs_ ? abs(v) : v.real());
				break;
			case NumberType::INT:
				*(int*)ptr_ = int(abs_ ? abs(v) : v.real());
				break;
			case NumberType::UINT:
				*(unsigned int*)ptr_ = (unsigned int)(abs_ ? abs(v) : v.real());
				break;
			case NumberType::SHORT:
				*(short*)ptr_ = short(abs_ ? abs(v) : v.real());
				break;
			case NumberType::USHORT:
				*(unsigned short*)ptr_ = (unsigned short)(abs_ ? abs(v) : v.real());
			}
			return *this;
		}
		template <typename T>
		NumRef& operator=(T v)
		{
			switch (type_) {
			case NumberType::CXDOUBLE:
				*(complex_double_t*)ptr_ = complex_double_t(v);
				break;
			case NumberType::CXFLOAT:
				*(complex_float_t*)ptr_ = complex_float_t(v);
				break;
			case NumberType::DOUBLE:
				*(double*)ptr_ = double(v);
				break;
			case NumberType::FLOAT:
				*(float*)ptr_ = float(v);
				break;
			case NumberType::INT:
				*(int*)ptr_ = int(v);
				break;
			case NumberType::UINT:
				*(unsigned int*)ptr_ = (unsigned int)v;
				break;
			case NumberType::SHORT:
				*(short*)ptr_ = short(v);
				break;
			case NumberType::USHORT:
				*(unsigned short*)ptr_ = (unsigned short)v;
			}
			return *this;
		}
        virtual NumberType::Type get_typeID() const
        {
            return NumberType::Type(type_);
        }

	private:
		bool abs_;
		int type_;
		void* ptr_;
	};

} // namespace sirf

#endif
