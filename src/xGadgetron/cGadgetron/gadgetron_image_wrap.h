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
\ingroup Gadgetron Image Wrapper
\brief Specification file for a wrapper class for ISMRMRD::Image.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef GADGETRON_IMAGE_WRAP_TYPE
#define GADGETRON_IMAGE_WRAP_TYPE

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include "cgadgetron_shared_ptr.h"
#include "xgadgetron_utilities.h"

#define IMAGE_PROCESSING_SWITCH(Type, Operation, Arguments, ...)\
	if (Type == ISMRMRD::ISMRMRD_USHORT)\
		Operation ((ISMRMRD::Image<unsigned short>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_SHORT)\
		Operation ((ISMRMRD::Image<short>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_UINT)\
		Operation ((ISMRMRD::Image<unsigned int>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_INT)\
		Operation ((ISMRMRD::Image<int>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_FLOAT)\
		Operation ((ISMRMRD::Image<float>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_DOUBLE)\
		Operation ((ISMRMRD::Image<double>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_CXFLOAT)\
		Operation ((ISMRMRD::Image< std::complex<float> >*) Arguments, \
			##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_CXDOUBLE)\
		Operation ((ISMRMRD::Image< std::complex<double> >*) Arguments, \
			##__VA_ARGS__);

#define IMAGE_PROCESSING_SWITCH_CONST(Type, Operation, Arguments, ...)\
	if (Type == ISMRMRD::ISMRMRD_USHORT)\
		Operation ((const ISMRMRD::Image<unsigned short>*) Arguments, \
			##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_SHORT)\
		Operation ((const ISMRMRD::Image<short>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_UINT)\
		Operation ((const ISMRMRD::Image<unsigned int>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_INT)\
		Operation ((const ISMRMRD::Image<int>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_FLOAT)\
		Operation ((const ISMRMRD::Image<float>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_DOUBLE)\
		Operation ((const ISMRMRD::Image<double>*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_CXFLOAT)\
		Operation ((const ISMRMRD::Image< std::complex<float> >*) Arguments, \
			##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_CXDOUBLE)\
		Operation ((const ISMRMRD::Image< std::complex<double> >*) Arguments, \
			##__VA_ARGS__);

typedef ISMRMRD::Image<complex_float_t> CFImage;
typedef ISMRMRD::Image<complex_double_t> CDImage;

namespace sirf {

	/**
\brief Wrapper for ISMRMRD::Image.

Eliminates the need for the image processing switch in the rest of the code.
*/
	class ImageWrap {
	public:
		ImageWrap(uint16_t type = 0, void* ptr_im = 0)
		{
			type_ = type;
			ptr_ = ptr_im;
		}
		ImageWrap(const ImageWrap& iw)
		{
			type_ = iw.type();
			IMAGE_PROCESSING_SWITCH(type_, copy_, iw.ptr_image());
		}
		~ImageWrap()
		{
			IMAGE_PROCESSING_SWITCH(type_, delete, ptr_);
		}
		int type() const
		{
			return type_;
		}
		void* ptr_image()
		{
			return ptr_;
		}
		const void* ptr_image() const
		{
			return ptr_;
		}
		size_t size() const
		{
			size_t s;
			IMAGE_PROCESSING_SWITCH_CONST(type_, get_size_, ptr_, s);
			return s;
		}
		ISMRMRD::ImageHeader& head()
		{
			IMAGE_PROCESSING_SWITCH(type_, return get_head_ref_, ptr_);
		}
		std::string attributes() const
		{
			std::string attr;
			IMAGE_PROCESSING_SWITCH_CONST(type_, get_attr_, ptr_, attr);
			return attr;
		}
		void set_imtype(ISMRMRD::ISMRMRD_ImageTypes imtype)
		{
			IMAGE_PROCESSING_SWITCH(type_, set_imtype_, ptr_, imtype);
		}
		size_t get_dim(int* dim) const
		{
			IMAGE_PROCESSING_SWITCH_CONST(type_, get_dim_, ptr_, dim);
			size_t n = dim[0];
			n *= dim[1];
			n *= dim[2];
			n *= dim[3];
			return n;
		}
		void get_data(float* data) const
		{
			IMAGE_PROCESSING_SWITCH_CONST(type_, get_data_, ptr_, data);
		}
		void set_data(const float* data) const
		{
			IMAGE_PROCESSING_SWITCH(type_, set_data_, ptr_, data);
		}
		void get_complex_data(complex_float_t* data) const
		{
			IMAGE_PROCESSING_SWITCH_CONST(type_, get_complex_data_, ptr_, data);
		}
		void set_complex_data(const complex_float_t* data) const
		{
			IMAGE_PROCESSING_SWITCH(type_, set_complex_data_, ptr_, data);
		}
		void write(ISMRMRD::Dataset& dataset) const
		{
			IMAGE_PROCESSING_SWITCH_CONST(type_, write_, ptr_, dataset);
		}
		void read(ISMRMRD::Dataset& dataset, const char* var, int ind)
		{
			IMAGE_PROCESSING_SWITCH(type_, read_, ptr_, dataset, var, ind, &ptr_);
		}
		void axpby(complex_float_t a, const ImageWrap& x, complex_float_t b)
		{
			IMAGE_PROCESSING_SWITCH(type_, axpby_, x.ptr_image(), a, b);
		}
		void multiply(const ImageWrap& x)
		{
			IMAGE_PROCESSING_SWITCH(type_, multiply_, x.ptr_image());
		}
		void divide(const ImageWrap& x)
		{
			IMAGE_PROCESSING_SWITCH(type_, divide_, x.ptr_image());
		}
		complex_float_t dot(const ImageWrap& iw) const
		{
			complex_float_t z;
			IMAGE_PROCESSING_SWITCH_CONST(type_, dot_, iw.ptr_image(), &z);
			return z;
		}
		float norm() const
		{
			float r;
			IMAGE_PROCESSING_SWITCH_CONST(type_, norm_, ptr_, &r);
			return r;
		}
		float diff(ImageWrap& iw) const
		{
			float s;
			IMAGE_PROCESSING_SWITCH_CONST(type_, diff_, iw.ptr_image(), &s);
			return s;
		}

		//void get_cmplx_data(float* re, float* im) const;

		//void set_cmplx_data(const float* re, const float* im) const;

	private:
		int type_;
		void* ptr_;

		ImageWrap& operator=(const ImageWrap& iw)
		{
			type_ = iw.type();
			IMAGE_PROCESSING_SWITCH(type_, copy_, iw.ptr_image());
			return *this;
		}

		template<typename T>
		void copy_(const ISMRMRD::Image<T>* ptr_im)
		{
			type_ = ptr_im->getDataType();
			ptr_ = (void*)new ISMRMRD::Image<T>(*ptr_im);
		}

		template<typename T>
		ISMRMRD::ImageHeader& get_head_ref_(ISMRMRD::Image<T>* ptr_im)
		{
			return ptr_im->getHead();
		}

		template<typename T>
		void set_imtype_(ISMRMRD::Image<T>* ptr_im, ISMRMRD::ISMRMRD_ImageTypes type)
		{
			ptr_im->setImageType(type);
		}

		template<typename T>
		void get_size_(const ISMRMRD::Image<T>* ptr_im, size_t& size) const
		{
			size = ptr_im->getDataSize();
		}

		template<typename T>
		void get_attr_(const ISMRMRD::Image<T>* ptr_im, std::string& attr) const
		{
			ptr_im->getAttributeString(attr);
		}

		template<typename T>
		void write_
			(const ISMRMRD::Image<T>* ptr_im, ISMRMRD::Dataset& dataset) const
		{
			//std::cout << "appending image..." << std::endl;
			const ISMRMRD::Image<T>& im = *ptr_im;
			std::stringstream ss;
			ss << "image_" << im.getHead().image_series_index;
			std::string image_varname = ss.str();
			{
				Mutex mtx;
				mtx.lock();
				dataset.appendImage(image_varname, im);
				mtx.unlock();
			}
		}

		template<typename T>
		void read_
			(const ISMRMRD::Image<T>* ptr,
			ISMRMRD::Dataset& dataset, const char* var, int index,
			void** ptr_ptr)
		{
			ISMRMRD::Image < T >* ptr_im = new ISMRMRD::Image < T > ;
			*ptr_ptr = (void*)ptr_im;
			ISMRMRD::Image<T>& im = *ptr_im;
			dataset.readImage(var, index, im);
			//int status = ismrmrd_read_image(&dataset, var, (uint32_t)index, &(im.im));
		}

		template<typename T>
		void get_dim_(const ISMRMRD::Image<T>* ptr_im, int* dim) const
		{
			const ISMRMRD::Image<T>& im = *(const ISMRMRD::Image<T>*)ptr_im;
			dim[0] = im.getMatrixSizeX();
			dim[1] = im.getMatrixSizeY();
			dim[2] = im.getMatrixSizeZ();
			dim[3] = im.getNumberOfChannels();
		}

		template<typename T>
		void get_data_(const ISMRMRD::Image<T>* ptr_im, float* data) const
		{
			const ISMRMRD::Image<T>& im = *ptr_im;
			const T* ptr = im.getDataPtr();
			size_t n = im.getNumberOfDataElements();
			for (size_t i = 0; i < n; i++)
				data[i] = std::real(ptr[i]);
		}

		template<typename T>
		void set_data_(ISMRMRD::Image<T>* ptr_im, const float* data) const
		{
			ISMRMRD::Image<T>& im = *ptr_im;
			T* ptr = im.getDataPtr();
			size_t n = im.getNumberOfDataElements();
			for (size_t i = 0; i < n; i++)
				ptr[i] = (T)data[i];
		}

		template<typename T>
		void get_complex_data_
			(const ISMRMRD::Image<T>* ptr_im, complex_float_t* data) const
		{
			const ISMRMRD::Image<T>& im = *ptr_im;
			const T* ptr = im.getDataPtr();
			size_t n = im.getNumberOfDataElements();
			for (size_t i = 0; i < n; i++)
				data[i] = complex_float_t(ptr[i]);
		}

		template<typename T>
		void set_complex_data_
			(ISMRMRD::Image<T>* ptr_im, const complex_float_t* data) const
		{
			ISMRMRD::Image<T>& im = *ptr_im;
			ISMRMRD::ISMRMRD_DataTypes type = im.getDataType();
			T* ptr = im.getDataPtr();
			size_t n = im.getNumberOfDataElements();
			if (type == ISMRMRD::ISMRMRD_CXFLOAT || type == ISMRMRD::ISMRMRD_CXDOUBLE)
				for (size_t i = 0; i < n; i++)
					xGadgetronUtilities::convert_complex(data[i], ptr[i]);
			//ptr[i] = data[i];
			else
				for (size_t i = 0; i < n; i++)
					xGadgetronUtilities::convert_complex(data[i], ptr[i]);
					//ptr[i] = std::real(data[i]);
		}

		template<typename T>
		void axpby_
			(const ISMRMRD::Image<T>* ptr_x, complex_float_t a, complex_float_t b)
		{
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)ptr_;
			const T* i;
			T* j;
			size_t ii = 0;
			size_t n = ptr_x->getNumberOfDataElements();
			if (b == complex_float_t(0.0))
				for (i = ptr_x->getDataPtr(), j = ptr_y->getDataPtr(); ii < n;
					i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				xGadgetronUtilities::convert_complex(a*u, *j);
			}
			else
				for (i = ptr_x->getDataPtr(), j = ptr_y->getDataPtr(); ii < n;
					i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				xGadgetronUtilities::convert_complex(a*u + b*v, *j);
			}
		}

		template<typename T>
		void multiply_(const ISMRMRD::Image<T>* ptr_x)
		{
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)ptr_;
			const T* i;
			T* j;
			size_t ii = 0;
			size_t n = ptr_x->getNumberOfDataElements();
			for (i = ptr_x->getDataPtr(), j = ptr_y->getDataPtr(); ii < n;
				i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				xGadgetronUtilities::convert_complex(u*v, *j);
			}
		}

		template<typename T>
		void divide_(const ISMRMRD::Image<T>* ptr_x)
		{
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)ptr_;
			const T* i;
			T* j;
			size_t ii = 0;
			size_t n = ptr_x->getNumberOfDataElements();
			for (i = ptr_x->getDataPtr(), j = ptr_y->getDataPtr(); ii < n;
				i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				// TODO: check for zero denominator
				xGadgetronUtilities::convert_complex(u / v, *j);
			}
		}

		template<typename T>
		void dot_(const ISMRMRD::Image<T>* ptr_im, complex_float_t *z) const
		{
			const ISMRMRD::Image<T>* ptr = (const ISMRMRD::Image<T>*)ptr_;
			const T* i;
			const T* j;
			*z = 0;
			size_t ii = 0;
			size_t n = ptr_im->getNumberOfDataElements();
			for (i = ptr->getDataPtr(), j = ptr_im->getDataPtr(); ii < n;
				i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				*z += std::conj(v) * u;
			}
		}

		template<typename T>
		void norm_(const ISMRMRD::Image<T>* ptr, float *r) const
		{
			const T* i;
			*r = 0;
			size_t ii = 0;
			size_t n = ptr->getNumberOfDataElements();
			for (i = ptr->getDataPtr(); ii < n; i++, ii++) {
				complex_float_t a = (complex_float_t)*i;
				*r += std::abs(std::conj(a) * a);
			}
			*r = std::sqrt(*r);
		}

		template<typename T>
		void diff_(const ISMRMRD::Image<T>* ptr_im, float *s) const
		{
			const ISMRMRD::Image<T>* ptr = (const ISMRMRD::Image<T>*)ptr_;
			const T* i;
			const T* j;
			*s = 0;
			size_t ii = 0;
			size_t n = ptr_im->getNumberOfDataElements();
			for (i = ptr->getDataPtr(), j = ptr_im->getDataPtr(); ii < n;
				i++, j++, ii++) {
				complex_float_t a = (complex_float_t)*i;
				complex_float_t b = (complex_float_t)*j;
				*s += (float)std::abs(b - a);
			}
		}
	};

}

#endif