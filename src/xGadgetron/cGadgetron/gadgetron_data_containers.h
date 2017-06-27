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

#ifndef GADGETRON_DATA_CONTAINERS
#define GADGETRON_DATA_CONTAINERS

#include <complex>
#include <fstream>
#include <map>
#include <thread>
#include <chrono>
#include <condition_variable>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include "ismrmrd_fftw.h"

#define TO_BE_IGNORED(acq) \
	(!(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) && \
	!(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING) && \
	!(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_MEASUREMENT) && \
	(acq).flags() >= (1 << (ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT - 1)))

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
		Operation ((ISMRMRD::Image< std::complex<float> >*) Arguments, ##__VA_ARGS__);\
	else if (Type == ISMRMRD::ISMRMRD_CXDOUBLE)\
		Operation ((ISMRMRD::Image< std::complex<double> >*) Arguments, ##__VA_ARGS__);

#define IMAGE_PROCESSING_SWITCH_CONST(Type, Operation, Arguments, ...)\
	if (Type == ISMRMRD::ISMRMRD_USHORT)\
		Operation ((const ISMRMRD::Image<unsigned short>*) Arguments, ##__VA_ARGS__);\
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
		Operation ((const ISMRMRD::Image< std::complex<float> >*) Arguments, ##__VA_ARGS__);\
		else if (Type == ISMRMRD::ISMRMRD_CXDOUBLE)\
		Operation ((const ISMRMRD::Image< std::complex<double> >*) Arguments, ##__VA_ARGS__);

typedef ISMRMRD::Image<complex_float_t> CFImage;
typedef ISMRMRD::Image<complex_double_t> CDImage;

class xGadgetronUtilities {
public:
	static long long milliseconds()
	{
		auto now = std::chrono::system_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
		return (long long)ms.count();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, unsigned short& t)
	{
		t = (unsigned short)z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, short& t)
	{
		t = (short)z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, unsigned int& t)
	{
		t = (unsigned int)z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, int& t)
	{
		t = (int)z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, float& t)
	{
		t = (float)z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, complex_float_t& t)
	{
		t = (complex_float_t)z;
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, double& t)
	{
		t = (double)z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, complex_double_t& t)
	{
		t = (complex_double_t)z;
	}

};

class Mutex {
public:
	Mutex() 
	{
		init_();
	}
	boost::mutex& operator()() 
	{
		return *sptr_mutex_.get();
	}
	boost::shared_ptr<boost::mutex> sptr() 
	{
		return sptr_mutex_;
	}
	void lock()
	{
		sptr_mutex_->lock();
	}
	void unlock()
	{
		sptr_mutex_->unlock();
	}
private:
	static boost::shared_ptr<boost::mutex> sptr_mutex_;
	static void init_() 
	{
		static bool initialized = false;
		if (!initialized) {
			sptr_mutex_ = boost::shared_ptr<boost::mutex>(new boost::mutex);
			initialized = true;
		}
	}
};

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
	ISMRMRD::ImageHeader* ptr_head()
	{
		ISMRMRD::ImageHeader** h;
		IMAGE_PROCESSING_SWITCH(type_, get_head_ptr_, ptr_, h);
		return *h;
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
	void get_data(double* data) const
	{
		IMAGE_PROCESSING_SWITCH_CONST(type_, get_data_, ptr_, data);
	}
	void write(ISMRMRD::Dataset& dataset) const
	{
		IMAGE_PROCESSING_SWITCH_CONST(type_, write_, ptr_, dataset);
	}
	void axpby(complex_double_t a, const ImageWrap& x, complex_double_t b)
	{
		IMAGE_PROCESSING_SWITCH(type_, axpby_, x.ptr_image(), a, b);
	}
	complex_double_t dot(const ImageWrap& iw) const
	{
		complex_double_t z;
		IMAGE_PROCESSING_SWITCH_CONST(type_, dot_, iw.ptr_image(), &z);
		return z;
	}
	double norm() const
	{
		double r;
		IMAGE_PROCESSING_SWITCH_CONST(type_, norm_, ptr_, &r);
		return r;
	}
	float diff(ImageWrap& iw) const
	{
		float s;
		IMAGE_PROCESSING_SWITCH_CONST(type_, diff_, iw.ptr_image(), &s);
		return s;
	}

	void get_cmplx_data(complex_float_t* data) const
	{
		const CFImage& im = *(const CFImage*)ptr_;
		long long int n = im.getMatrixSizeX();
		n *= im.getMatrixSizeY();
		n *= im.getMatrixSizeZ();
		n *= im.getNumberOfChannels();
		const complex_float_t* ptr = im.getDataPtr();
		for (long long int i = 0; i < n; i++)
			data[i] = ptr[i];
	}

	void get_cmplx_data(double* re, double* im) const
	{
		int dim[4];
		size_t n = get_dim(dim);
		if (type_ == ISMRMRD::ISMRMRD_CXFLOAT) {
			const CFImage& img = *(const CFImage*)ptr_;
			const complex_float_t* ptr = img.getDataPtr();
			for (size_t i = 0; i < n; i++) {
				complex_float_t z = ptr[i];
				re[i] = std::real(z);
				im[i] = std::imag(z);
			}
		}
		else if (type_ == ISMRMRD::ISMRMRD_CXDOUBLE) {
			const CDImage& img = *(const CDImage*)ptr_;
			const complex_double_t* ptr = img.getDataPtr();
			for (size_t i = 0; i < n; i++) {
				complex_double_t z = ptr[i];
				re[i] = std::real(z);
				im[i] = std::imag(z);
			}
		}
		else {
			get_data(re);
			for (size_t i = 0; i < n; i++)
				im[i] = 0;
		}
	}

	void set_cmplx_data(const double* re, const double* im) const
	{
		int dim[4];
		size_t n = get_dim(dim);
		if (type_ == ISMRMRD::ISMRMRD_CXFLOAT) {
			CFImage& img = *(CFImage*)ptr_;
			complex_float_t* ptr = img.getDataPtr();
			for (size_t i = 0; i < n; i++)
				ptr[i] = std::complex<float>((float)re[i], (float)im[i]);
		}
		else if (type_ == ISMRMRD::ISMRMRD_CXDOUBLE) {
			CDImage& img = *(CDImage*)ptr_;
			complex_double_t* ptr = img.getDataPtr();
			for (size_t i = 0; i < n; i++)
				ptr[i] = std::complex<double>(re[i], im[i]);
		}
	}

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
	void get_head_ptr_(ISMRMRD::Image<T>* ptr_im, ISMRMRD::ImageHeader** h)
	{
		*h = &(ptr_im->getHead());
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
	void get_dim_(const ISMRMRD::Image<T>* ptr_im, int* dim) const
	{
		const ISMRMRD::Image<T>& im = *(const ISMRMRD::Image<T>*)ptr_im;
		dim[0] = im.getMatrixSizeX();
		dim[1] = im.getMatrixSizeY();
		dim[2] = im.getMatrixSizeZ();
		dim[3] = im.getNumberOfChannels();
	}

	template<typename T>
	void get_data_(const ISMRMRD::Image<T>* ptr_im, double* data) const
	{
		const ISMRMRD::Image<T>& im = *ptr_im;
		const T* ptr = im.getDataPtr();
		size_t n = im.getNumberOfDataElements();
		for (size_t i = 0; i < n; i++)
			data[i] = std::real(ptr[i]);
	}

	template<typename T>
	void axpby_
		(const ISMRMRD::Image<T>* ptr_x, complex_double_t a, complex_double_t b)
	{
		ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)ptr_;
		const T* i;
		T* j;
		size_t ii = 0;
		size_t n = ptr_x->getNumberOfDataElements();
		if (b == complex_double_t(0.0))
			for (i = ptr_x->getDataPtr(), j = ptr_y->getDataPtr(); ii < n;
				i++, j++, ii++) {
			complex_double_t u = (complex_double_t)*i;
			xGadgetronUtilities::convert_complex(a*u, *j);
		}
		else
			for (i = ptr_x->getDataPtr(), j = ptr_y->getDataPtr(); ii < n;
				i++, j++, ii++) {
			complex_double_t u = (complex_double_t)*i;
			complex_double_t v = (complex_double_t)*j;
			xGadgetronUtilities::convert_complex(a*u + b*v, *j);
		}
	}

	template<typename T>
	void dot_(const ISMRMRD::Image<T>* ptr_im, complex_double_t *z) const
	{
		const ISMRMRD::Image<T>* ptr = (const ISMRMRD::Image<T>*)ptr_;
		const T* i;
		const T* j;
		*z = 0;
		size_t ii = 0;
		size_t n = ptr_im->getNumberOfDataElements();
		for (i = ptr->getDataPtr(), j = ptr_im->getDataPtr(); ii < n;
			i++, j++, ii++) {
			complex_double_t u = (complex_double_t)*i;
			complex_double_t v = (complex_double_t)*j;
			*z += std::conj(v) * u;
		}
	}

	template<typename T>
	void norm_(const ISMRMRD::Image<T>* ptr, double *r) const
	{
		const T* i;
		*r = 0;
		size_t ii = 0;
		size_t n = ptr->getNumberOfDataElements();
		for (i = ptr->getDataPtr(); ii < n; i++, ii++) {
			complex_double_t a = (complex_double_t)*i;
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
			complex_double_t a = (complex_double_t)*i;
			complex_double_t b = (complex_double_t)*j;
			*s += (float)std::abs(b - a);
		}
	}
};

class aDataContainer {
public:
	virtual ~aDataContainer() {}
	virtual boost::shared_ptr<aDataContainer> new_data_container() = 0;
	virtual unsigned int items() = 0;
	virtual double norm() = 0;
	virtual complex_double_t dot(aDataContainer& dc) = 0;
	virtual void axpby(
		complex_double_t a, const aDataContainer& a_x,
		complex_double_t b, const aDataContainer& a_y) = 0;
};

class AcquisitionsContainer : public aDataContainer {
public:
	AcquisitionsContainer() : ordered_(false), index_(0) {}
	virtual ~AcquisitionsContainer()
	{
		if (index_)
			delete[] index_;
	}

	std::string parameters() const
	{
		return par_;
	}
	void set_parameters(std::string par)
	{
		par_ = par;
	}
	bool undersampled() const
	{
		ISMRMRD::IsmrmrdHeader header;
		ISMRMRD::deserialize(par_.c_str(), header);
		ISMRMRD::Encoding e = header.encoding[0];
		return e.parallelImaging.is_present() &&
			e.parallelImaging().accelerationFactor.kspace_encoding_step_1 > 1;
	}
	int get_acquisitions_dimensions(size_t ptr_dim)
	{
		ISMRMRD::Acquisition acq;
		int* dim = (int*)ptr_dim;

		int na = number();
		int ms, ns;
		int mc, nc;
		int my, ny;
		int slice = 0;
		int y = 0;
		// assume all dimensions (samples, coils [, acqs per slice]) regular
		int nrd = ordered() ? 3 : 2;
		// number of regular readouts
		int nrr = 0;
		//int not_reg = 0;
		for (; y < na;){
			for (; y < na && ordered();){
				get_acquisition(y, acq);
				if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
					break;
				y++;
			}
			if (y >= na)
				break;
			ny = 0;
			for (; y < na; y++) {
				get_acquisition(y, acq);
				if (TO_BE_IGNORED(acq)) // not a regular acquisition
					continue;
				ns = acq.number_of_samples();
				nc = acq.active_channels();
				nrr += ns*nc;
				if (slice == 0) {
					ms = ns;
					mc = nc;
				}
				else {
					if (ms != ns)
						nrd = 0;
					else if (mc != nc && nrd > 1)
						nrd = 1;
				}
				ny++;
				if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE) && ordered())
					break;
			}
			if (slice == 0) {
				my = ny;
			}
			else {
				if (my != ny && nrd > 2) {
					nrd = 2;
				}
				//if (my != ny)
				//	not_reg = 1;
			}
			slice++;
		}

		int reg_size = 1;
		if (nrd > 0) {
			dim[0] = ms;
			reg_size *= ms;
		}
		if (nrd > 1) {
			dim[1] = mc;
			reg_size *= mc;
		}
		if (nrd > 2) {
			dim[2] = my;
			reg_size *= my;
		}
		dim[nrd] = nrr / reg_size;
		return nrd;
		//dim[0] = acq.number_of_samples();
		//dim[1] = acq.active_channels();
		//dim[2] = my; // e.reconSpace.matrixSize.y;
		//dim[3] = slice;
		//return not_reg;
	}

	int set_acquisitions_data
		 (boost::shared_ptr<AcquisitionsContainer> sptr_ac, 
		 int na, int nc, int ns, const double* re, const double* im)
	{
		sptr_ac->set_parameters(par_);
		sptr_ac->write_parameters();
		sptr_ac->ordered_ = ordered();
		ISMRMRD::Acquisition acq;
		int ma = number();
		for (int a = 0, i = 0; a < ma; a++) {
			get_acquisition(a, acq);
			if (TO_BE_IGNORED(acq) && ma > na) {
				std::cout << "ignoring acquisition " << a << '\n';
				continue;
			}
			unsigned int mc = acq.active_channels();
			unsigned int ms = acq.number_of_samples();
			if (mc != nc || ms != ns)
				return -1;
			//for (size_t c = 0; c < nc; c++) {
			//	for (size_t s = 0; s < ns; s++, i++) {
			//		acq.data(s, c) = complex_float_t((float)re[i], (float)im[i]);
			//	}
			//}
			for (int c = 0; c < nc; c++)
				for (int s = 0; s < ns; s++, i++)
					acq.data(s, c) = complex_float_t((float)re[i], (float)im[i]);
			sptr_ac->append_acquisition(acq);
		}
		return 0;
	}

	void get_acquisitions_flags(unsigned int n, int* flags)
	{
		ISMRMRD::Acquisition acq;
		unsigned int na = number();
		for (unsigned int a = 0, i = 0; a < na; a++) {
			get_acquisition(a, acq);
			if (TO_BE_IGNORED(acq) && n < na) {
				std::cout << "ignoring acquisition " << a << '\n';
				continue;
			}
			flags[i++] = (int)acq.flags();
		}
	}

	unsigned int get_acquisitions_data(unsigned int slice, double* re, double* im)
	{
		ISMRMRD::Acquisition acq;
		unsigned int na = number();
		unsigned int n = 0;
		if (slice >= na) {
			for (unsigned int a = 0, i = 0; a < na; a++) {
				get_acquisition(a, acq);
				if (TO_BE_IGNORED(acq) && slice > na) {
					std::cout << "ignoring acquisition " << a << '\n';
					continue;
				}
				n++;
				unsigned int nc = acq.active_channels();
				unsigned int ns = acq.number_of_samples();
				for (unsigned int c = 0; c < nc; c++) {
					for (unsigned int s = 0; s < ns; s++, i++) {
						complex_float_t z = acq.data(s, c);
						re[i] = std::real(z);
						im[i] = std::imag(z);
					}
				}
			}
			return n;
		}
		int* dim = new int[3];
		size_t ptr_dim = (size_t)dim;
		get_acquisitions_dimensions(ptr_dim);
		unsigned int ny = dim[2]; //e.reconSpace.matrixSize.y;
		//unsigned int ny = dim[1]; //e.reconSpace.matrixSize.y;
		delete[] dim;
		unsigned int y = 0;
		for (; y + ny*slice < na;){
			get_acquisition(y + ny*slice, acq);
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
				break;
			y++;
		}
		for (; y + ny*slice < na; n++) {
			get_acquisition(y + ny*slice, acq);
			unsigned int nc = acq.active_channels();
			unsigned int ns = acq.number_of_samples();
			for (unsigned int c = 0; c < nc; c++) {
				for (unsigned int s = 0; s < ns; s++) {
					complex_float_t z = acq.data(s, c);
					re[s + ns*(n + ny*c)] = std::real(z);
					im[s + ns*(n + ny*c)] = std::imag(z);
				}
			}
			y++;
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
				break;
		}
		return n;
	}

	static void axpby
		(complex_double_t a, const ISMRMRD::Acquisition& acq_x, 
		complex_double_t b, ISMRMRD::Acquisition& acq_y)
	{
		complex_float_t* px;
		complex_float_t* py;
		for (px = acq_x.data_begin(), py = acq_y.data_begin();
			px != acq_x.data_end() && py != acq_y.data_end(); px++, py++) {
			if (b == complex_double_t(0.0))
				*py = a*complex_double_t(*px);
			else 
				*py = a*complex_double_t(*px) + b*complex_double_t(*py);
		}
	}
	static complex_double_t dot
		(const ISMRMRD::Acquisition& acq_a, const ISMRMRD::Acquisition& acq_b)
	{
		complex_float_t* pa;
		complex_float_t* pb;
		complex_double_t z = 0;
		for (pa = acq_a.data_begin(), pb = acq_b.data_begin();
			pa != acq_a.data_end() && pb != acq_b.data_end(); pa++, pb++) {
			z += std::conj(*pb) * (*pa);
		}
		return z;
	}
	static double norm(const ISMRMRD::Acquisition& acq_a)
	{
		complex_float_t* pa;
		double r = 0;
		for (pa = acq_a.data_begin(); pa != acq_a.data_end(); pa++) {
			complex_float_t z = std::conj(*pa) * (*pa);
			r += z.real();
		}
		r = sqrt(r);
		return r;
	}
	static float diff
		(const ISMRMRD::Acquisition& acq_a, const ISMRMRD::Acquisition& acq_b)
	{
		complex_float_t* pa;
		complex_float_t* pb;
		complex_float_t z = 0;
		float s = 0.0;
		float sa = 0.0;
		float sb = 0.0;
		for (pa = acq_a.data_begin(), pb = acq_b.data_begin();
			pa != acq_a.data_end() && pb != acq_b.data_end(); pa++, pb++) {
			float ta = std::abs(*pa);
			float tb = std::abs(*pb);
			sa += pow(std::abs(*pa), 2);
			sb += pow(std::abs(*pb), 2);
			z += std::conj(*pb) * (*pa);
		}
		z /= sb;
		sa = std::sqrt(sa);
		sb = std::sqrt(sb);
		for (pa = acq_a.data_begin(), pb = acq_b.data_begin();
			pa != acq_a.data_end() && pb != acq_b.data_end(); pa++, pb++) {
			s += pow(std::abs(*pa - *pb * z), 2);
		}
		s = std::sqrt(s);
		s /= sa;
		return s;
	}

	virtual unsigned int number() = 0;
	virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
	virtual void append_acquisition(ISMRMRD::Acquisition& acq) = 0;
	virtual void copy_parameters(const AcquisitionsContainer& ac) = 0;
	virtual void write_parameters() = 0;
	virtual 
		boost::shared_ptr<AcquisitionsContainer> new_acquisitions_container() = 0;

	virtual void axpby(
		complex_double_t a, const aDataContainer& a_x,
		complex_double_t b, const aDataContainer& a_y)
	{
		AcquisitionsContainer& x = (AcquisitionsContainer&)a_x;
		AcquisitionsContainer& y = (AcquisitionsContainer&)a_y;
		int m = x.number();
		int n = y.number();
		ISMRMRD::Acquisition ax;
		ISMRMRD::Acquisition ay;
		for (int i = 0, j = 0; i < n && j < m;) {
			y.get_acquisition(i, ay);
			x.get_acquisition(j, ax);
			if (TO_BE_IGNORED(ay)) {
				std::cout << i << " ignored (ay)\n";
				i++;
				continue;
			}
			if (TO_BE_IGNORED(ax)) {
				std::cout << j << " ignored (ax)\n";
				j++;
				continue;
			}
			AcquisitionsContainer::axpby(a, ax, b, ay);
			append_acquisition(ay);
			i++;
			j++;
		}
	}
	virtual complex_double_t dot(aDataContainer& dc)
	{
		AcquisitionsContainer& other = (AcquisitionsContainer&)dc;
		int n = number();
		int m = other.number();
		complex_double_t z = 0;
		ISMRMRD::Acquisition a;
		ISMRMRD::Acquisition b;
		for (int i = 0, j = 0; i < n && j < m;) {
			get_acquisition(i, a);
			if (TO_BE_IGNORED(a)) {
				i++;
				continue;
			}
			other.get_acquisition(j, b);
			if (TO_BE_IGNORED(b)) {
				j++;
				continue;
			}
			z += AcquisitionsContainer::dot(a, b);
			i++;
			j++;
		}
		return z;
	}
	virtual double norm()
	{
		int n = number();
		double r = 0;
		ISMRMRD::Acquisition a;
		for (int i = 0; i < n; i++) {
			get_acquisition(i, a);
			if (TO_BE_IGNORED(a)) {
				continue;
			}
			double s = AcquisitionsContainer::norm(a);
			r += s*s;
		}
		return sqrt(r);
	}

	float diff(AcquisitionsContainer& other)
	{
		int n = number();
		int m = other.number();
		float smax = 0.0;
		float save = 0.0;
		ISMRMRD::Acquisition a;
		ISMRMRD::Acquisition b;
		for (int i = 0; i < n && i < m; i++) {
			get_acquisition(i, a);
			other.get_acquisition(i, b);
			float s = AcquisitionsContainer::diff(a, b);
			smax = std::max(smax, s);
			save += s*s;
		}
		save = sqrt(save / std::min(n, m));
		return save;
	}

	void order()
	{
		class triple {
		public:
			triple(int i, int f, int s, int t) : 
				ind(i), first(f), second(s), third(t) {}
			triple(const triple& t) : 
				ind(t.ind), first(t.first), second(t.second), third(t.third) {}
			int ind;
			int first;
			int second;
			int third;
		};

		int na = number();
		std::vector<triple> t;
		ISMRMRD::Acquisition acq;
		for (int i = 0; i < na; i++) {
			get_acquisition(i, acq);
			int rep = acq.idx().repetition;
			int slice = acq.idx().slice;
			int phase = acq.idx().kspace_encode_step_1;
			t.push_back(triple(i, rep, slice, phase));
		}
		std::stable_sort(t.begin(), t.end(), 
			[](triple a, triple b) { return b.first > a.first; });

		int i = 0;
		int j = i;
		while (i < na) {
			std::vector<triple> ts;
			for (; j < na; j++) {
				if (t[j].first != t[i].first)
					break;
			}
			for (int k = i; k < j; k++)
				ts.push_back(triple(t[k]));
			std::stable_sort(ts.begin(), ts.end(),
				[](triple a, triple b) { return b.second > a.second; });
			for (int k = i; k < j; k++)
				t[k] = ts[k - i];
			i = j;
		}

		if (index_)
			delete[] index_;
		index_ = new int[na];

		i = 0;
		j = i;
		while (i < na) {
			std::vector<triple> ts;
			for (; j < na; j++) {
				if (t[j].first != t[i].first || t[j].second != t[i].second)
					break;
			}
			for (int k = i; k < j; k++)
				ts.push_back(triple(t[k]));
			std::sort(ts.begin(), ts.end(),
				[](triple a, triple b) { return b.third > a.third; });
			for (int k = i; k < j; k++)
				index_[k] = ts[k - i].ind;
			i = j;
		}
		ordered_ = true;
	}
	bool ordered() const
	{
		return ordered_;
		//return (bool)index_;
	}
	int index(int i)
	{
		if (index_ && i >= 0 && i < (int)number())
			return index_[i];
		else
			return i;
	}

protected:
	bool ordered_;
	int* index_;
	std::string par_;
};

class AcquisitionsFile : public AcquisitionsContainer {
public:
	AcquisitionsFile
		(std::string filename, bool create_file = false, bool own_file = false)
	{
		own_file_ = own_file;
		filename_ = filename;
		int ndim = 0;
		Mutex mtx;
		mtx.lock();
		dataset_ = boost::shared_ptr<ISMRMRD::Dataset>
			(new ISMRMRD::Dataset(filename.c_str(), "/dataset", create_file));
		if (!create_file) {
			dataset_->readHeader(par_);
		}
		mtx.unlock();
	}
	~AcquisitionsFile() {
		dataset_.reset();
		if (own_file_) {
			Mutex mtx;
			mtx.lock();
			std::remove(filename_.c_str());
			mtx.unlock();
		}

	}
	virtual unsigned int items()
	{
		Mutex mtx;
		mtx.lock();
		unsigned int na = dataset_->getNumberOfAcquisitions();
		mtx.unlock();
		return na;
	}
	virtual unsigned int number()
	{
		return items();
	}
	virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
	{
		int ind = index(num);
		Mutex mtx;
		mtx.lock();
		dataset_->readAcquisition(ind, acq);
		//dataset_->readAcquisition(index(num), acq); // ??? does not work!
		mtx.unlock();
	}
	virtual void append_acquisition(ISMRMRD::Acquisition& acq)
	{
		Mutex mtx;
		mtx.lock();
		dataset_->appendAcquisition(acq);
		mtx.unlock();
	}
	virtual void copy_parameters(const AcquisitionsContainer& ac) 
	{
		par_ = ac.parameters();
		Mutex mtx;
		mtx.lock();
		dataset_->writeHeader(par_);
		mtx.unlock();
	}
	virtual void write_parameters()
	{
		Mutex mtx;
		mtx.lock();
		dataset_->writeHeader(par_);
		mtx.unlock();
	}
	virtual boost::shared_ptr<aDataContainer> new_data_container()
	{
		AcquisitionsFile* ptr_ac = acqs_scratch_file_(filename_);
		ptr_ac->set_parameters(par_);
		ptr_ac->write_parameters();
		boost::shared_ptr<aDataContainer> sptr_ac(ptr_ac);
		return sptr_ac;
	}
	virtual boost::shared_ptr<AcquisitionsContainer> new_acquisitions_container()
	{
		boost::shared_ptr<AcquisitionsContainer> 
			sptr_ac(acqs_scratch_file_(filename_));
		sptr_ac->set_parameters(par_);
		sptr_ac->write_parameters();
		return sptr_ac;
	}

private:
	bool own_file_;
	std::string filename_;
	boost::shared_ptr<ISMRMRD::Dataset> dataset_;

	AcquisitionsFile* acqs_scratch_file_(std::string name)
	{
		static int calls = 0;
		char buff[32];
		long long int ms = xGadgetronUtilities::milliseconds();
		calls++;
		sprintf(buff, "_%d_%lld.h5", calls, ms);
		boost::replace_all(name, ".h5", buff);
		//std::cout << "new acquisitions file: " << name << std::endl;
		return new AcquisitionsFile(name, true, true);
	}
};

class AcquisitionsList : public AcquisitionsContainer {
public:
	virtual unsigned int number()
	{
		return (int)acqs_.size();
	}
private:
	std::list<boost::shared_ptr<ISMRMRD::Acquisition> > acqs_;
};

class ImagesContainer : public aDataContainer {
public:
	virtual unsigned int number() = 0;
	virtual int types() = 0;
	virtual void count(int i) = 0;
	virtual boost::shared_ptr<ImageWrap> sptr_image_wrap(unsigned int im_num) = 0;
	virtual boost::shared_ptr<const ImageWrap> sptr_image_wrap
		(unsigned int im_num) const = 0;
	virtual ImageWrap& image_wrap(unsigned int im_num) = 0;
	virtual const ImageWrap& image_wrap(unsigned int im_num) const = 0;
	virtual void append(int image_data_type, void* ptr_image) = 0;
	virtual void append(const ImageWrap& iw) = 0;
	virtual void get_image_dimensions(unsigned int im_num, int* dim) = 0;
	virtual void get_image_data_as_double_array
		(unsigned int im_num, double* data) = 0;
	virtual void get_image_data_as_complex_array
		(unsigned int im_num, complex_float_t* data) = 0;
	virtual void get_images_data_as_double_array(double* data) const = 0;
	virtual void get_images_data_as_complex_array
		(double* re, double* im) const = 0;
	virtual void set_complex_images_data(const double* re, const double* im) = 0;
	virtual void write(std::string filename, std::string groupname) = 0;
	virtual boost::shared_ptr<ImagesContainer> new_images_container() = 0;
	virtual boost::shared_ptr<ImagesContainer>
		clone(unsigned int inc = 1, unsigned int off = 0) = 0;
	virtual boost::shared_ptr<ImagesContainer>
		clone(const char* attr, const char* target) = 0;

	virtual int image_data_type(unsigned int im_num) const
	{
		return image_wrap(im_num).type();
	}

	virtual void axpby(
		complex_double_t a, const aDataContainer& a_x,
		complex_double_t b, const aDataContainer& a_y)
	{
		ImagesContainer& x = (ImagesContainer&)a_x;
		ImagesContainer& y = (ImagesContainer&)a_y;
		ImageWrap w(x.image_wrap(0));
		complex_double_t zero(0.0, 0.0);
		complex_double_t one(1.0, 0.0);
		for (unsigned int i = 0; i < x.number() && i < y.number(); i++) {
			const ImageWrap& u = x.image_wrap(i);
			const ImageWrap& v = y.image_wrap(i);
			w.axpby(a, u, zero);
			w.axpby(b, v, one);
			append(w);
		}
	}
	virtual complex_double_t dot(aDataContainer& dc)
	{
		ImagesContainer& ic = (ImagesContainer&)dc;
		complex_double_t z = 0;
		for (unsigned int i = 0; i < number() && i < ic.number(); i++) {
			const ImageWrap& u = image_wrap(i);
			const ImageWrap& v = ic.image_wrap(i);
			z += u.dot(v);
		}
		return z;
	}
	virtual double norm()
	{
		double r = 0;
		for (unsigned int i = 0; i < number(); i++) {
			const ImageWrap& u = image_wrap(i);
			double s = u.norm();
			r += s*s;
		}
		r = sqrt(r);
		return r;
	}

	void get_image_data_as_cmplx_array
		(unsigned int im_num, double* re, double* im)
	{
		ImageWrap& iw = image_wrap(im_num);
		iw.get_cmplx_data(re, im);
	}

	void set_image_to_real_conversion(int type)
	{
		for (unsigned int i = 0; i < number(); i++) {
			ImageWrap& u = image_wrap(i);
			u.set_imtype((ISMRMRD::ISMRMRD_ImageTypes)type);
		}
	}
};

class ImagesList : public ImagesContainer {
public:
	ImagesList() : images_(), nimages_(0)
	{
	}
	ImagesList(const ImagesList& list, const char* attr, const char* target)
	{
#ifdef _MSC_VER
		std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#endif
		for (i = list.images_.begin(); i != list.images_.end(); i++) {
			const boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			std::string atts = sptr_iw->attributes();
			ISMRMRD::MetaContainer mc;
			ISMRMRD::deserialize(atts.c_str(), mc);
			std::string value = mc.as_str(attr);
			if (boost::iequals(value, target))
				append(*sptr_iw);
		}
	}
	ImagesList(const ImagesList& list, unsigned int inc = 1, unsigned int off = 0)
	{
		int n = 0;
		unsigned int j = 0;
#ifdef _MSC_VER
		std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#endif
		for (i = list.images_.begin(); i != list.images_.end() && j < off; i++, j++);

		for (; i != list.images_.end(); i++, j++) {
			if ((j - off) % inc)
				continue;
			const boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			append(*sptr_iw);
			n++;
		}
		nimages_ = n;
	}
	virtual unsigned int items()
	{
		return (unsigned int)images_.size();
	}
	virtual unsigned int number()
	{
		return (unsigned int)images_.size();
	}
	virtual int types()
	{
		if (nimages_ > 0)
			return (int)(images_.size() / nimages_);
		else
			return 1;
	}
	virtual void count(int i)
	{
		if (i > nimages_)
			nimages_ = i;
	}
	virtual void append(int image_data_type, void* ptr_image)
	{
		images_.push_back(boost::shared_ptr<ImageWrap>
			(new ImageWrap(image_data_type, ptr_image)));
	}
	virtual void append(const ImageWrap& iw)
	{
		images_.push_back(boost::shared_ptr<ImageWrap>(new ImageWrap(iw)));
	}
	virtual boost::shared_ptr<ImageWrap> sptr_image_wrap(unsigned int im_num)
	{
#ifdef _MSC_VER
		std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#endif
		unsigned int count = 0;
		for (i = images_.begin(); 
			i != images_.end() && count < im_num && count < images_.size() - 1; i++)
			count++;
		return *i;
	}
	virtual boost::shared_ptr<const ImageWrap> sptr_image_wrap
		(unsigned int im_num) const
	{
#ifdef _MSC_VER
		std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#endif
		unsigned int count = 0;
		for (i = images_.begin(); i != images_.end() && count < im_num; i++)
			count++;
		return *i;
	}
	virtual ImageWrap& image_wrap(unsigned int im_num)
	{
		boost::shared_ptr<ImageWrap> sptr_iw = sptr_image_wrap(im_num);
		return *sptr_iw;
	}
	virtual const ImageWrap& image_wrap(unsigned int im_num) const
	{
		const boost::shared_ptr<const ImageWrap>& sptr_iw = sptr_image_wrap(im_num);
		return *sptr_iw;
	}
	virtual void write(std::string filename, std::string groupname)
	{
		if (images_.size() < 1)
			return;
		Mutex mtx;
		mtx.lock();
		ISMRMRD::Dataset dataset(filename.c_str(), groupname.c_str());
		mtx.unlock();
#ifdef _MSC_VER
		std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#endif
		for (i = images_.begin(); i != images_.end(); i++) {
			boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			ImageWrap& iw = *sptr_iw;
			iw.write(dataset);
		}
	}
	virtual void get_image_dimensions(unsigned int im_num, int* dim)
	{
		if (im_num >= images_.size())
			dim[0] = dim[1] = dim[2] = dim[3] = 0;
		ImageWrap& iw = image_wrap(im_num);
		iw.get_dim(dim);
		//std::string attr = iw.attributes();
		//ISMRMRD::MetaContainer mc;
		//ISMRMRD::deserialize(attr.c_str(), mc);
		//std::cout << mc.as_str("GADGETRON_DataRole") << '\n';
		//std::cout << attr << '\n';
	}
	virtual void get_image_data_as_double_array(unsigned int im_num, double* data)
	{
		ImageWrap& iw = image_wrap(im_num);
		iw.get_data(data);
	}
	virtual void get_image_data_as_complex_array
		(unsigned int im_num, complex_float_t* data)
	{
		ImageWrap& iw = image_wrap(im_num);
		iw.get_cmplx_data(data);
	}
	virtual void get_images_data_as_double_array(double* data) const
	{
#ifdef _MSC_VER
		std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#endif
		int dim[4];
		for (i = images_.begin(); i != images_.end(); i++) {
			const boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			ImageWrap& iw = *sptr_iw;
			iw.get_data(data);
			iw.get_dim(dim);
			size_t size = dim[0];
			size *= dim[1];
			size *= dim[2];
			size *= dim[3];
			data += size;
		}
	}
	virtual void get_images_data_as_complex_array(double* re, double* im) const
	{
#ifdef _MSC_VER
		std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#endif
		int dim[4];
		for (i = images_.begin(); i != images_.end(); i++) {
			const boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			ImageWrap& iw = *sptr_iw;
			iw.get_dim(dim);
			size_t size = dim[0];
			size *= dim[1];
			size *= dim[2];
			size *= dim[3];
			int type = iw.type();
			if (type == ISMRMRD::ISMRMRD_CXFLOAT || type == ISMRMRD::ISMRMRD_CXDOUBLE)
				iw.get_cmplx_data(re, im);
			else {
				iw.get_data(re);
				for (int i = 0; i < size; i++)
					im[i] = 0;
			}
			re += size;
			im += size;
		}
	}
	virtual void set_complex_images_data(const double* re, const double* im)
	{
#ifdef _MSC_VER
		std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#endif
		int dim[4];
		for (i = images_.begin(); i != images_.end(); i++) {
			boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			ImageWrap& iw = *sptr_iw;
			size_t n = iw.get_dim(dim);
			iw.set_cmplx_data(re, im);
			re += n;
			im += n;
		}
	}
	virtual boost::shared_ptr<aDataContainer> new_data_container()
	{
		return boost::shared_ptr<aDataContainer>((aDataContainer*)new ImagesList());
	}
	virtual boost::shared_ptr<ImagesContainer> new_images_container()
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesList());
	}
	virtual boost::shared_ptr<ImagesContainer>
		clone(const char* attr, const char* target)
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesList(*this, attr, target));
	}
	virtual boost::shared_ptr<ImagesContainer>
		clone(unsigned int inc = 1, unsigned int off = 0)
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesList(*this, inc, off));
	}

private:
	std::list<boost::shared_ptr<ImageWrap> > images_;
	int nimages_;
};

class CoilData {
public:
	virtual ~CoilData() {}
	virtual void get_dim(int* dim) const = 0;
	virtual void get_data(double* re, double* im) const = 0;
	virtual void set_data(const double* re, const double* im) = 0;
	virtual void get_data(complex_float_t* data) const = 0;
	virtual void set_data(const complex_float_t* data) = 0;
	virtual void get_data_abs(double* v) const = 0;
	virtual complex_float_t& operator()(int x, int y, int z, int c) = 0;
};

class CoilDataAsCFImage : public CoilData {
public:
	CoilDataAsCFImage(uint16_t nx = 0, uint16_t ny = 1, uint16_t nz = 1, uint16_t nc = 1) :
		img_(nx, ny, nz, nc)
	{
	}
	virtual complex_float_t& operator()(int x, int y, int z, int c)
	{
		return img_(x, y, z, c);
	}
	ISMRMRD::Image < complex_float_t >& image()
	{
		return img_;
	}
	const ISMRMRD::Image < complex_float_t >& image() const
	{
		return img_;
	}
	virtual void get_dim(int* dim) const
	{
		dim[0] = img_.getMatrixSizeX();
		dim[1] = img_.getMatrixSizeY();
		dim[2] = img_.getMatrixSizeZ();
		dim[3] = img_.getNumberOfChannels();
	}
	virtual void get_data(double* re, double* im) const
	{
		size_t n = img_.getNumberOfDataElements();
		const complex_float_t* ptr = img_.getDataPtr();
		for (size_t i = 0; i < n; i++) {
			complex_float_t z = ptr[i];
			re[i] = std::real(z);
			im[i] = std::imag(z);
		}
	}
	virtual void set_data(const double* re, const double* im)
	{
		size_t n = img_.getNumberOfDataElements();
		complex_float_t* ptr = img_.getDataPtr();
		for (size_t i = 0; i < n; i++)
			ptr[i] = complex_float_t((float)re[i], (float)im[i]);
	}
	virtual void get_data(complex_float_t* data) const
	{
		memcpy(data, img_.getDataPtr(), img_.getDataSize());
	}
	virtual void set_data(const complex_float_t* data)
	{
		memcpy(img_.getDataPtr(), data, img_.getDataSize());
	}
	virtual void get_data_abs(double* v) const
	{
		size_t n = img_.getNumberOfDataElements();
		const complex_float_t* ptr = img_.getDataPtr();
		for (size_t i = 0; i < n; i++) {
			complex_float_t z = ptr[i];
			v[i] = std::abs(z);
		}
	}
private:
	ISMRMRD::Image < complex_float_t > img_;
};

class CoilDataContainer : public aDataContainer {
public:
	virtual double norm()
	{
		return 0.0;
	}
	virtual complex_double_t dot(aDataContainer& dc)
	{
		return complex_double_t(0.0, 0.0);
	}
	virtual void axpby(
		complex_double_t a, const aDataContainer& a_x,
		complex_double_t b, const aDataContainer& a_y)
	{
		return;
	}
	void get_dim(int slice, int* dim) //const
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.get_dim(dim);
	}
	void get_data(int slice, double* re, double* im) //const
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.get_data(re, im);
	}
	void set_data(int slice, double* re, double* im)
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.set_data(re, im);
	}
	void get_data(int slice, complex_float_t* data) //const
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.get_data(data);
	}
	void set_data(int slice, complex_float_t* data)
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.set_data(data);
	}
	void get_data_abs(int slice, double* v) //const
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.get_data_abs(v);
	}
	virtual void append(boost::shared_ptr<CoilData> sptr_csm) = 0;
	virtual CoilData& operator()(int slice) = 0;
	//virtual const CoilData& operator()(int slice) const = 0;
};

class CoilDataList {
public:
	unsigned int items()
	{
		return (unsigned int)list_.size();
	}
	CoilData& data(int slice)
	{
#ifdef _MSC_VER
		std::list<boost::shared_ptr<CoilData> >::const_iterator i;
#else
		typename std::list<boost::shared_ptr<CoilData> >::const_iterator i;
#endif
		int count = 0;
		for (i = list_.begin();
			i != list_.end() && count < slice && count < (int)list_.size() - 1; i++)
			count++;
		return **i;
	}
	virtual void append(boost::shared_ptr<CoilData> sptr_cd)
	{
		list_.push_back(sptr_cd);
	}
private:
	std::list< boost::shared_ptr<CoilData> > list_;
};

class CoilImagesContainer : public CoilDataContainer {
public:
	virtual CoilData& operator()(int slice) = 0;
	virtual void compute(AcquisitionsContainer& ac)
	{
		std::string par;
		ISMRMRD::IsmrmrdHeader header;
		ISMRMRD::Acquisition acq;
		par = ac.parameters();
		ISMRMRD::deserialize(par.c_str(), header);
		ac.get_acquisition(0, acq);
		encoding_ = header.encoding[0];

		ISMRMRD::Encoding e = header.encoding[0];
		bool parallel = e.parallelImaging.is_present() &&
			e.parallelImaging().accelerationFactor.kspace_encoding_step_1 > 1;
		unsigned int nx = e.reconSpace.matrixSize.x;
		unsigned int ny = e.reconSpace.matrixSize.y;
		unsigned int nc = acq.active_channels();
		unsigned int readout = acq.number_of_samples();

		//std::cout << nx << ' ' << ny << ' ' << nc << ' ' << readout << '\n';
		//if (e.parallelImaging.is_present()) {
		//	std::cout << "parallel imaging present\n";
		//	std::cout << "acceleration factors: " 
		//		<< e.parallelImaging().accelerationFactor.kspace_encoding_step_1 << ' '
		//		<< e.parallelImaging().accelerationFactor.kspace_encoding_step_2 << '\n';
		//}
		//else
		//	std::cout << "parallel imaging not present\n";

		int nmap = 0;
		std::cout << "map ";

		for (unsigned int na = 0; na < ac.number();) {

			std::cout << ++nmap << ' ' << std::flush;

			std::vector<size_t> ci_dims;
			ci_dims.push_back(readout);
			ci_dims.push_back(ny);
			ci_dims.push_back(nc);
			ISMRMRD::NDArray<complex_float_t> ci(ci_dims);
			memset(ci.getDataPtr(), 0, ci.getDataSize());

			int y = 0;
			for (;;){
				ac.get_acquisition(na + y, acq);
				if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
					break;
				y++;
			}
			for (;;) {
				ac.get_acquisition(na + y, acq);
				int yy = acq.idx().kspace_encode_step_1;
				//if (!e.parallelImaging.is_present() ||
				if ( !parallel ||
					acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) ||
					acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING)) {
					for (unsigned int c = 0; c < nc; c++) {
						for (unsigned int s = 0; s < readout; s++) {
							ci(s, yy, c) = acq.data(s, c);
						}
					}
				}
				y++;
				if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
					break;
			}
			na += y;

			ifft2c(ci);

			boost::shared_ptr<CoilData>
				sptr_ci(new CoilDataAsCFImage(readout, ny, 1, nc));
			CFImage& coil_im = (*(CoilDataAsCFImage*)sptr_ci.get()).image();
			memcpy(coil_im.getDataPtr(), ci.getDataPtr(), ci.getDataSize());
			append(sptr_ci);
		}
		std::cout << '\n';
	}
	ISMRMRD::Encoding encoding() const
	{
		return encoding_;
	}
protected:
	ISMRMRD::Encoding encoding_;
};

class CoilImagesList : public CoilImagesContainer, public CoilDataList {
public:
	virtual boost::shared_ptr<aDataContainer> new_data_container()
	{
		return boost::shared_ptr<aDataContainer>
			((aDataContainer*)new CoilImagesList());
	}
	virtual unsigned int items()
	{
		return CoilDataList::items();
	}
	virtual CoilData& operator()(int slice)
	{
		return data(slice);
	}
	virtual void append(boost::shared_ptr<CoilData> sptr_cd)
	{
		CoilDataList::append(sptr_cd);
	}
};

template<class CoilDataType>
class CoilSensitivitiesContainerTemplate : public CoilDataContainer {
public:
	void set_csm_smoothness(int s)
	{
		csm_smoothness_ = s;
	}
	virtual CoilData& operator()(int slice) = 0;

	virtual void compute(AcquisitionsContainer& ac)
	{
		//if (!ac.ordered())
		//	ac.order();

		CoilImagesList cis;
		cis.compute(ac);
		compute(cis);
	}

	virtual void compute(CoilImagesContainer& cis)
	{

		ISMRMRD::Encoding e = cis.encoding();
		unsigned int nx = e.reconSpace.matrixSize.x;
		unsigned int ny = e.reconSpace.matrixSize.y;
		int dim[4];
		cis(0).get_dim(dim);
		unsigned int readout = dim[0];
		unsigned int nc = dim[3];

		std::vector<size_t> cm_dims;
		cm_dims.push_back(readout);
		cm_dims.push_back(ny);
		cm_dims.push_back(nc);
		ISMRMRD::NDArray<complex_float_t> cm(cm_dims);

		std::vector<size_t> csm_dims;
		csm_dims.push_back(nx);
		csm_dims.push_back(ny);
		csm_dims.push_back(1);
		csm_dims.push_back(nc);
		ISMRMRD::NDArray<complex_float_t> csm(csm_dims);

		std::vector<size_t> img_dims;
		img_dims.push_back(nx);
		img_dims.push_back(ny);
		ISMRMRD::NDArray<float> img(img_dims);

		unsigned int nmap = 0;

		std::cout << "map ";
		for (nmap = 1; nmap <= cis.items(); nmap++) {
			std::cout << nmap << ' ' << std::flush;
			cis(nmap - 1).get_data(cm.getDataPtr());
			CoilData* ptr_img = new CoilDataType(nx, ny, 1, nc);
			boost::shared_ptr<CoilData> sptr_img(ptr_img);
			compute_(cm, img, csm);
			ptr_img->set_data(csm.getDataPtr());
			append(sptr_img);
		}
		std::cout << '\n';
	}

	void append_csm
		(int nx, int ny, int nz, int nc, const double* re, const double* im)
	{
		CoilData* ptr_img = new CoilDataType(nx, ny, nz, nc);
		boost::shared_ptr<CoilData> sptr_img(ptr_img);
		ptr_img->set_data(re, im);
		append(sptr_img);
	}

protected:
	int csm_smoothness_;

private:
	void compute_(
		ISMRMRD::NDArray<complex_float_t>& cm,
		ISMRMRD::NDArray<float>& img,
		ISMRMRD::NDArray<complex_float_t>& csm
		)
	{
		int ndims = cm.getNDim();
		const size_t* dims = cm.getDims();
		unsigned int readout = (unsigned int)dims[0];
		unsigned int ny = (unsigned int)dims[1];
		unsigned int nc = (unsigned int)dims[2];
		unsigned int nx = (unsigned int)img.getDims()[0];

		std::vector<size_t> cm0_dims;
		cm0_dims.push_back(nx);
		cm0_dims.push_back(ny);
		cm0_dims.push_back(nc);

		ISMRMRD::NDArray<complex_float_t> cm0(cm0_dims);
		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++) {
					uint16_t xout = x + (readout - nx) / 2;
					cm0(x, y, c) = cm(xout, y, c);
				}
			}
		}

		int* object_mask = new int[nx*ny*nc];
		memset(object_mask, 0, nx*ny*nc*sizeof(int));

		ISMRMRD::NDArray<complex_float_t> w(cm0);

		float* ptr_img = img.getDataPtr();
		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				double r = 0.0;
				for (unsigned int c = 0; c < nc; c++) {
					float s = std::abs(cm0(x, y, c));
					r += s*s;
				}
				img(x, y) = (float)std::sqrt(r);
			}
		}

		float noise = max_(5, 5, ptr_img) + (float)1e-6*max_(nx, ny, ptr_img);
		mask_noise_(nx, ny, ptr_img, noise, object_mask);
		cleanup_mask_(nx, ny, object_mask, 0, 2, 0);
		cleanup_mask_(nx, ny, object_mask, 0, 3, 0);
		cleanup_mask_(nx, ny, object_mask, 0, 4, 0);

		for (int i = 0; i < csm_smoothness_; i++)
			smoothen_(nx, ny, nc, cm0.getDataPtr(), w.getDataPtr(), object_mask);

		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				double r = 0.0;
				for (unsigned int c = 0; c < nc; c++) {
					float s = std::abs(cm0(x, y, c));
					r += s*s;
				}
				img(x, y) = (float)std::sqrt(r);
			}
		}

		for (unsigned int y = 0, i = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++, i++) {
				double r = img(x, y);
				float s;
				if (r != 0.0)
					s = (float)(1.0 / r);
				else
					s = 0.0;
				complex_float_t z(s, 0.0);
				for (unsigned int c = 0; c < nc; c++) {
					csm(x, y, 0, c) = cm0(x, y, c) * z;
				}
			}
		}

		delete[] object_mask;

	}

	float max_(int nx, int ny, float* u)
	{
		float r = 0.0;
		int i = 0;
		for (int iy = 0; iy < ny; iy++)
			for (int ix = 0; ix < nx; ix++, i++) {
				float t = fabs(u[i]);
				if (t > r)
					r = t;
			}
		return r;
	}
	void mask_noise_
		(int nx, int ny, float* u, float noise, int* mask)
	{
		int i = 0;
		for (int iy = 0; iy < ny; iy++)
			for (int ix = 0; ix < nx; ix++, i++) {
				float t = fabs(u[i]);
				mask[i] = (t > noise);
			}
	}
	void cleanup_mask_(int nx, int ny, int* mask, int bg, int minsz, int ex)
	{
		int ll, il;
		int* listx = new int[nx*ny];
		int* listy = new int[nx*ny];
		int* inlist = new int[nx*ny];
		std::memset(inlist, 0, nx*ny*sizeof(int));
		for (int iy = 0, i = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++, i++) {
				if (mask[i] == bg)
					continue;
				bool skip = false;
				ll = 1;
				listx[0] = ix;
				listy[0] = iy;
				inlist[i] = 1;
				il = 0;
				while (il < ll && ll < minsz) {
					int lx = listx[il];
					int ly = listy[il];
					int l = ll + ex;
					for (int jy = -l; jy <= l; jy++) {
						for (int jx = -l; jx <= l; jx++) {
							int kx = lx + jx;
							int ky = ly + jy;
							if (kx < 0 || kx >= nx)
								continue;
							if (ky < 0 || ky >= ny)
								continue;
							int j = kx + ky*nx;
							if (inlist[j])
								continue;
							if (mask[j] != bg) {
								listx[ll] = kx;
								listy[ll] = ky;
								inlist[j] = 1;
								ll++;
							}
						}
					}
					il++;
				}
				if (il == ll)
					mask[i] = bg;
				for (il = 0; il < ll; il++) {
					int lx = listx[il];
					int ly = listy[il];
					int j = lx + ly*nx;
					inlist[j] = 0;
				}
			}
		}
		delete[] listx;
		delete[] listy;
		delete[] inlist;
	}
	void smoothen_
		(int nx, int ny, int nz, 
		complex_float_t* u, complex_float_t* v, 
		int* obj_mask)
	{
		const complex_float_t ONE(1.0, 0.0);
		const complex_float_t TWO(2.0, 0.0);
		for (int iz = 0, i = 0; iz < nz; iz++)
			for (int iy = 0, k = 0; iy < ny; iy++)
				for (int ix = 0; ix < nx; ix++, i++, k++) {
					//if (edge_mask[k]) {
					//	v[i] = u[i];
					//	continue;
					//}
					int n = 0;
					complex_float_t r(0.0, 0.0);
					complex_float_t s(0.0, 0.0);
					for (int jy = -1; jy <= 1; jy++)
						for (int jx = -1; jx <= 1; jx++) {
							if (ix + jx < 0 || ix + jx >= nx)
								continue;
							if (iy + jy < 0 || iy + jy >= ny)
								continue;
							int j = i + jx + jy*nx;
							int l = k + jx + jy*nx;
							if (i != j && obj_mask[l]) { // && !edge_mask[l]) {
								n++;
								r += ONE;
								s += u[j];
							}
						}
					if (n > 0)
						v[i] = (u[i] + s / r) / TWO;
					else
						v[i] = u[i];
				}
		memcpy(u, v, nx*ny*nz*sizeof(complex_float_t));
	}
};

typedef CoilSensitivitiesContainerTemplate<CoilDataAsCFImage> CoilSensitivitiesContainer;

class CoilSensitivitiesAsImages : public CoilSensitivitiesContainer, 
	public CoilDataList {
public:
	CoilSensitivitiesAsImages()
	{
		csm_smoothness_ = 0;
	}
	CoilSensitivitiesAsImages(const char* file)
	{
		Mutex mtx;
		mtx.lock();
		ISMRMRD::Dataset csm_file(file, "dataset");
		int nm = csm_file.getNumberOfImages("csm");
		mtx.unlock();
		for (int i = 0; i < nm; i++) {
			boost::shared_ptr<CoilData> sptr_img(new CoilDataAsCFImage);
			mtx.lock();
			CFImage& csm = (*(CoilDataAsCFImage*)sptr_img.get()).image();
			csm_file.readImage("csm", i, csm);
			mtx.unlock();
			append(sptr_img);
		}
		csm_smoothness_ = 0;
	}

	virtual boost::shared_ptr<aDataContainer> new_data_container()
	{
		return boost::shared_ptr<aDataContainer>
			((aDataContainer*)new CoilSensitivitiesAsImages());
	}

	virtual unsigned int items()
	{
		return CoilDataList::items();
	}
	virtual CoilData& operator()(int slice)
	{
		return data(slice);
	}
	virtual void append(boost::shared_ptr<CoilData> sptr_cd)
	{
		CoilDataList::append(sptr_cd);
	}

};

#endif
