/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
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

/*!
\file
\ingroup Gadgetron Image Wrapper
\brief Specification file for a wrapper class for ISMRMRD::Image.

\author Evgueni Ovtchinnikov
\author SyneRBI
*/

#ifndef GADGETRON_IMAGE_WRAP_TYPE
#define GADGETRON_IMAGE_WRAP_TYPE

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include "sirf/common/ANumRef.h"
#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/xgadgetron_utilities.h"

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
		class Iterator {
		public:
			Iterator(int type, void* data, unsigned int dsize, size_t n) :
				type_(type), ptr_((char*)data), dsize_(dsize), n_(n), i_(0),
				ref_(data, type)
			{}
			Iterator(const Iterator& iter)
			{
				type_ = iter.type_;
				ptr_ = iter.ptr_;
				dsize_ = iter.dsize_;
				n_ = iter.n_;
				i_ = iter.i_;
				sptr_iter_ = iter.sptr_iter_;
				ref_.copy(iter.ref_);
			}
			bool operator!=(const Iterator& i) const
			{
				return ptr_ != i.ptr_;
			}
			bool operator==(const Iterator& i) const
			{
				return ptr_ == i.ptr_;
			}
			Iterator& operator++()
			{
				if (i_ >= n_)
					throw std::out_of_range("cannot advance out-of-range iterator");
				i_++;
				ptr_ += dsize_;
				return *this;
			}
			// too inefficient
			//virtual Iterator& operator++(int)
			//{
			//	//sptr_iter_.reset(new Iterator(type_, ptr_, dsize_, n_));
			//	sptr_iter_.reset(new Iterator(*this));
			//	if (i_ >= n_)
			//		throw std::out_of_range("cannot advance out-of-range iterator");
			//	i_++;
			//	ptr_ += dsize_;
			//	return *sptr_iter_;
			//}
			NumRef& operator*()
			{
				if (i_ >= n_)
					throw std::out_of_range
					("cannot dereference out-of-range iterator");
				ref_.set_ptr(ptr_);
				return ref_;
			}
			Iterator& operator=(const Iterator& iter)
			{
				type_ = iter.type_;
				ptr_ = iter.ptr_;
				dsize_ = iter.dsize_;
				n_ = iter.n_;
				i_ = iter.i_;
				sptr_iter_ = iter.sptr_iter_;
				ref_.copy(iter.ref_);
				return *this;
			}
		private:
			int type_;
			char* ptr_;
			unsigned int dsize_;
			size_t n_;
			size_t i_;
			gadgetron::shared_ptr<Iterator> sptr_iter_;
			NumRef ref_;
		};
		class Iterator_const {
		public:
			Iterator_const(int type, void* data, unsigned int dsize, size_t n) :
				type_(type), ptr_((char*)data), dsize_(dsize), n_(n), i_(0),
				ref_(data, type)
			{}
			Iterator_const(const Iterator_const& iter)
			{
				type_ = iter.type_;
				ptr_ = iter.ptr_;
				dsize_ = iter.dsize_;
				n_ = iter.n_;
				i_ = iter.i_;
				sptr_iter_ = iter.sptr_iter_;
				ref_.copy(iter.ref_);
			}
			bool operator!=(const Iterator_const& i) const
			{
				return ptr_ != i.ptr_;
			}
			bool operator==(const Iterator_const& i) const
			{
				return ptr_ == i.ptr_;
			}
			Iterator_const& operator++()
			{
				if (i_ >= n_)
					throw std::out_of_range("cannot advance out-of-range iterator");
				i_++;
				ptr_ += dsize_;
				return *this;
			}
			//virtual Iterator_const& operator++(int)
			//{
			//	//sptr_iter_.reset(new Iterator_const(type_, ptr_, dsize_, n_));
			//	sptr_iter_.reset(new Iterator_const(*this));
			//	if (i_ >= n_)
			//		throw std::out_of_range("cannot advance out-of-range iterator");
			//	i_++;
			//	ptr_ += dsize_;
			//	return *sptr_iter_;
			//}
			const NumRef& operator*() const
			{
				if (i_ >= n_) {
					std::cout << i_ << ' ' << n_ << '\n';
					throw std::out_of_range
						("cannot dereference out-of-range iterator");
				}
				ref_.set_ptr(ptr_);
				return ref_;
			}
			Iterator_const& operator=(const Iterator_const& iter)
			{
				type_ = iter.type_;
				ptr_ = iter.ptr_;
				dsize_ = iter.dsize_;
				n_ = iter.n_;
				i_ = iter.i_;
				sptr_iter_ = iter.sptr_iter_;
				ref_.copy(iter.ref_);
				return *this;
			}
		private:
			int type_;
			char* ptr_;
			unsigned int dsize_;
			size_t n_;
			size_t i_;
			gadgetron::shared_ptr<Iterator_const> sptr_iter_;
			mutable NumRef ref_;
		};
		ImageWrap(uint16_t type, void* ptr_im)
		{
			type_ = type;
			ptr_ = ptr_im;
		}
		ImageWrap(uint16_t type, ISMRMRD::Dataset& dataset, const char* var, int index)
		{
			type_ = type;
			IMAGE_PROCESSING_SWITCH(type_, read_, ptr_, dataset, var, index, &ptr_);
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
		Iterator& begin()
		{
			size_t n;
			unsigned int dsize;
			char* ptr;
			IMAGE_PROCESSING_SWITCH
			(type_, get_data_parameters_, ptr_, &n, &dsize, &ptr);
			begin_.reset(new Iterator(type_, ptr, dsize, n));
			return *begin_;
		}
		Iterator_const& begin_const() const
		{
			size_t n;
			unsigned int dsize;
			char* ptr;
			IMAGE_PROCESSING_SWITCH_CONST
			(type_, get_data_parameters_, ptr_, &n, &dsize, &ptr);
			begin_const_.reset(new Iterator_const(type_, ptr, dsize, n));
			return *begin_const_;
		}
		Iterator& end()
		{
			size_t n;
			unsigned int dsize;
			char* ptr;
			IMAGE_PROCESSING_SWITCH
			(type_, get_data_parameters_, ptr_, &n, &dsize, &ptr);
			end_.reset(new Iterator(type_, ptr + n*dsize, dsize, n));
			return *end_;
		}
		Iterator_const& end_const() const
		{
			size_t n;
			unsigned int dsize;
			char* ptr;
			IMAGE_PROCESSING_SWITCH_CONST
			(type_, get_data_parameters_, ptr_, &n, &dsize, &ptr);
			end_const_.reset(new Iterator_const(type_, ptr + n*dsize, dsize, n));
			//std::cout << type_ << ' ' << n << ' ' << dsize << '\n';
			return *end_const_;
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
			//std::cout << "in get_data\n";
			//std::cout << "trying new image wrap iterator...\n";
			ImageWrap::Iterator_const i = begin_const();
			ImageWrap::Iterator_const stop = end_const();
			for (; i != stop; ++data, ++i) {
				*data = *i;
			}
			//IMAGE_PROCESSING_SWITCH_CONST(type_, get_data_, ptr_, data);
		}
		void set_data(const float* data)
		{
			//std::cout << "in set_data\n";
			for (ImageWrap::Iterator i = begin(); i != end(); ++i, ++data)
				*i = *data;
			//IMAGE_PROCESSING_SWITCH(type_, set_data_, ptr_, data);
		}
		void get_complex_data(complex_float_t* data) const
		{
			//std::cout << "in get_complex_data\n";
			//std::cout << "trying new const image wrap iterator...\n";
			ImageWrap::Iterator_const i = begin_const();
			ImageWrap::Iterator_const stop = end_const();
			for (; i != stop; ++data, ++i) {
				*data = (*i).complex_float();
			}
			//IMAGE_PROCESSING_SWITCH_CONST(type_, get_complex_data_, ptr_, data);
		}

        void set_complex_data(const complex_float_t* data)
		{
			//std::cout << "in set_complex_data\n";
			//std::cout << "trying new image wrap iterator...\n";
			ImageWrap::Iterator i = begin();
			ImageWrap::Iterator stop = end();
			for (; i != stop; ++i, ++data) {
				*i = *data;
			}
			//IMAGE_PROCESSING_SWITCH(type_, set_complex_data_, ptr_, data);
		}

        /// Get data type
        ISMRMRD::ISMRMRD_DataTypes get_data_type() const
        {
            ISMRMRD::ISMRMRD_DataTypes data_type;
            IMAGE_PROCESSING_SWITCH_CONST(type_, get_data_type_, ptr_, &data_type);
            return data_type;
        }
        /// Is the image wrap complex?
        bool is_complex() const
        {
            ISMRMRD::ISMRMRD_DataTypes data_type = this->get_data_type();
            if (data_type == ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT ||
                    data_type == ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXDOUBLE) {
                return true;
            }
            else
                return false;
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
		void axpby(complex_float_t a, const ImageWrap& x, complex_float_t b, 
			const ImageWrap& y)
		{
			IMAGE_PROCESSING_SWITCH(type_, axpby__, x.ptr_image(), a, y.ptr_image(), b);
		}
		void multiply(const ImageWrap& x)
		{
			IMAGE_PROCESSING_SWITCH(type_, multiply_, x.ptr_image());
		}
		void multiply(const ImageWrap& x, const ImageWrap& y)
		{
			IMAGE_PROCESSING_SWITCH(type_, multiply__, x.ptr_image(), y.ptr_image());
		}
		void divide(const ImageWrap& x)
		{
			IMAGE_PROCESSING_SWITCH(type_, divide_, x.ptr_image());
		}
		void divide(const ImageWrap& x, const ImageWrap& y)
		{
			IMAGE_PROCESSING_SWITCH(type_, divide__, x.ptr_image(), y.ptr_image());
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

	private:
		int type_;
		void* ptr_;
		mutable gadgetron::shared_ptr<Iterator> begin_;
		mutable gadgetron::shared_ptr<Iterator> end_;
		mutable gadgetron::shared_ptr<Iterator_const> begin_const_;
		mutable gadgetron::shared_ptr<Iterator_const> end_const_;

		ImageWrap& operator=(const ImageWrap& iw)
		{
			//type_ = iw.type();
			//IMAGE_PROCESSING_SWITCH(type_, copy_, iw.ptr_image());
			return *this;
		}

		template<typename T>
		void get_data_parameters_
			(const ISMRMRD::Image<T>* ptr_im, size_t *n, unsigned int *dsize,
			char** data) const
		{
			const ISMRMRD::Image<T>& im = *(const ISMRMRD::Image<T>*)ptr_im;
			unsigned int dim[4];
			dim[0] = im.getMatrixSizeX();
			dim[1] = im.getMatrixSizeY();
			dim[2] = im.getMatrixSizeZ();
			dim[3] = im.getNumberOfChannels();
			*data = (char*)im.getDataPtr();
			*dsize = sizeof(T);
			*n = dim[0];
			*n *= dim[1];
			*n *= dim[2];
			*n *= dim[3];
			//*n = ptr_im->getDataSize()/(*dsize);
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
		void get_data_type_(const ISMRMRD::Image<T>* ptr_im, ISMRMRD::ISMRMRD_DataTypes* data_type_ptr) const
		{
			const ISMRMRD::Image<T>& im = *static_cast<const ISMRMRD::Image<T>*>(ptr_im);
			*data_type_ptr = im.getDataType();
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
		void set_data_(ISMRMRD::Image<T>* ptr_im, const float* data)
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
			(ISMRMRD::Image<T>* ptr_im, const complex_float_t* data)
		{
			ISMRMRD::Image<T>& im = *ptr_im;
			ISMRMRD::ISMRMRD_DataTypes type = im.getDataType();
			T* ptr = im.getDataPtr();
			size_t n = im.getNumberOfDataElements();
			if (type == ISMRMRD::ISMRMRD_CXFLOAT || type == ISMRMRD::ISMRMRD_CXDOUBLE)
				for (size_t i = 0; i < n; i++)
					xGadgetronUtilities::convert_complex(data[i], ptr[i]);
			else
				for (size_t i = 0; i < n; i++)
					xGadgetronUtilities::convert_complex(data[i], ptr[i]);
		}

		template<typename T>
		void axpby_
			(const ISMRMRD::Image<T>* ptr_x, complex_float_t a, 
			complex_float_t b)
		{
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)ptr_;
			size_t nx = ptr_x->getNumberOfDataElements();
			size_t ny = ptr_y->getNumberOfDataElements();
			if (nx != ny)
				THROW("sizes mismatch in ImageWrap multiply");
			const T* i = ptr_x->getDataPtr();
			T* j = ptr_y->getDataPtr();
			if (b == complex_float_t(0.0))
				for (size_t ii = 0; ii < nx; i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				xGadgetronUtilities::convert_complex(a*u, *j);
			}
			else
				for (size_t ii = 0; ii < nx; i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				xGadgetronUtilities::convert_complex(a*u + b*v, *j);
			}
		}

		template<typename T>
		void axpby__
			(const ISMRMRD::Image<T>* ptr_x, complex_float_t a, 
			const void* vptr_y, complex_float_t b)
		{
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)vptr_y;
			ISMRMRD::Image<T>* ptr = (ISMRMRD::Image<T>*)ptr_;
			size_t n = ptr->getNumberOfDataElements();
			size_t nx = ptr_x->getNumberOfDataElements();
			size_t ny = ptr_y->getNumberOfDataElements();
			if (!(n == nx && n == ny))
				THROW("sizes mismatch in ImageWrap multiply");
			const T* i = ptr_x->getDataPtr();
			const T* j = ptr_y->getDataPtr();
			T* k = ptr->getDataPtr();
			if (b == complex_float_t(0.0))
				for (size_t ii = 0; ii < n; i++, k++, ii++) {
					complex_float_t u = (complex_float_t)*i;
					xGadgetronUtilities::convert_complex(a*u, *k);
			}
			else
				for (size_t ii = 0; ii < n; i++, j++, k++, ii++) {
					complex_float_t u = (complex_float_t)*i;
					complex_float_t v = (complex_float_t)*j;
					complex_float_t w = a*u + b*v;
					xGadgetronUtilities::convert_complex(w, *k);
			}
		}

		template<typename T>
		void multiply_(const ISMRMRD::Image<T>* ptr_x)
		{
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)ptr_;
			size_t nx = ptr_x->getNumberOfDataElements();
			size_t ny = ptr_y->getNumberOfDataElements();
			if (nx != ny)
				THROW("sizes mismatch in ImageWrap multiply");
			const T* i = ptr_x->getDataPtr();
			T* j = ptr_y->getDataPtr();
			size_t ii = 0;
			for (; ii < nx; i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				xGadgetronUtilities::convert_complex(u*v, *j);
			}
		}

		template<typename T>
		void multiply__(const ISMRMRD::Image<T>* ptr_x, const void* vptr_y)
		{
			ISMRMRD::Image<T>* ptr = (ISMRMRD::Image<T>*)ptr_;
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)vptr_y;
			size_t nx = ptr_x->getNumberOfDataElements();
			size_t ny = ptr_y->getNumberOfDataElements();
			size_t n = ptr->getNumberOfDataElements();
			if (!(n == nx && n == ny))
				THROW("sizes mismatch in ImageWrap multiply");
			const T* i = ptr_x->getDataPtr();
			const T* j = ptr_y->getDataPtr();
			T* k = ptr->getDataPtr();
			size_t ii = 0;
			for (; ii < n; i++, j++, k++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				xGadgetronUtilities::convert_complex(u*v, *k);
			}
		}

		template<typename T>
		void divide_(const ISMRMRD::Image<T>* ptr_x)
		{
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)ptr_;
			size_t nx = ptr_x->getNumberOfDataElements();
			size_t ny = ptr_y->getNumberOfDataElements();
			if (nx != ny)
				THROW("sizes mismatch in ImageWrap multiply");
			const T* i = ptr_x->getDataPtr();
			T* j = ptr_y->getDataPtr();
			size_t ii = 0;
			for (; ii < nx; i++, j++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				// TODO: check for zero denominator
				xGadgetronUtilities::convert_complex(v / u, *j);
			}
		}

		template<typename T>
		void divide__(const ISMRMRD::Image<T>* ptr_x, const void* vptr_y)
		{
			ISMRMRD::Image<T>* ptr = (ISMRMRD::Image<T>*)ptr_;
			ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)vptr_y;
			size_t nx = ptr_x->getNumberOfDataElements();
			size_t ny = ptr_y->getNumberOfDataElements();
			size_t n = ptr->getNumberOfDataElements();
			if (!(n == nx && n == ny))
				THROW("sizes mismatch in ImageWrap multiply");
			const T* i = ptr_x->getDataPtr();
			const T* j = ptr_y->getDataPtr();
			T* k = ptr->getDataPtr();
			size_t ii = 0;
			for (; ii < n; i++, j++, k++, ii++) {
				complex_float_t u = (complex_float_t)*i;
				complex_float_t v = (complex_float_t)*j;
				// TODO: check for zero denominator
				xGadgetronUtilities::convert_complex(u / v, *k);
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
