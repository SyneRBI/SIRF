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
\ingroup Gadgetron Data Containers
\brief Specification file for data container classes for Gadgetron data.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef GADGETRON_DATA_CONTAINERS
#define GADGETRON_DATA_CONTAINERS

#include <algorithm> // stable_sort
#include <array> // array
#include <numeric> // iota
#include <vector> // vector

namespace Multisort {

	template<typename T, size_t N>
	bool less(const std::array<T, N>& a, const std::array<T, N>& b)
	{
		for (unsigned int i = 0; i < N; i++) {
			if (a[i] < b[i])
				return true;
			if (a[i] > b[i])
				return false;
		}
		return false; // all equal
	}

	template<typename T, size_t N>
	void sort(std::vector<std::array<T, N> > v, int* index)
	{
		int n = v.size();
		std::iota(index, index + n, 0);
		std::stable_sort
			(index, index + n, [&v](int i, int j){return less(v[i], v[j]); });
	}

} // namespace Multisort

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
	void get_data(float* data) const
	{
		IMAGE_PROCESSING_SWITCH_CONST(type_, get_data_, ptr_, data);
	}
	void write(ISMRMRD::Dataset& dataset) const
	{
		IMAGE_PROCESSING_SWITCH_CONST(type_, write_, ptr_, dataset);
	}
	void axpby(complex_float_t a, const ImageWrap& x, complex_float_t b)
	{
		IMAGE_PROCESSING_SWITCH(type_, axpby_, x.ptr_image(), a, b);
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

	void get_cmplx_data(float* re, float* im) const;

	void set_cmplx_data(const float* re, const float* im) const;

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
	void get_data_(const ISMRMRD::Image<T>* ptr_im, float* data) const
	{
		const ISMRMRD::Image<T>& im = *ptr_im;
		const T* ptr = im.getDataPtr();
		size_t n = im.getNumberOfDataElements();
		for (size_t i = 0; i < n; i++)
			data[i] = std::real(ptr[i]);
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

template <typename T>
class aDataContainer {
public:
	virtual ~aDataContainer() {}
	virtual boost::shared_ptr<aDataContainer<T> > new_data_container() = 0;
	virtual unsigned int items() = 0;
	virtual float norm() = 0;
	virtual T dot(aDataContainer<T>& dc) = 0;
	virtual void axpby(
		T a, const aDataContainer<T>& a_x,
		T b, const aDataContainer<T>& a_y) = 0;
};

class AcquisitionsContainer : public aDataContainer<complex_float_t> {
public:
	AcquisitionsContainer() : ordered_(false), index_(0) {}
	virtual ~AcquisitionsContainer()
	{
		if (index_)
			delete[] index_;
	}

	static void axpby
	(complex_float_t a, const ISMRMRD::Acquisition& acq_x,
		complex_float_t b, ISMRMRD::Acquisition& acq_y);
	static complex_float_t dot
	(const ISMRMRD::Acquisition& acq_a, const ISMRMRD::Acquisition& acq_b);
	static float norm(const ISMRMRD::Acquisition& acq_a);
	static float diff
	(const ISMRMRD::Acquisition& acq_a, const ISMRMRD::Acquisition& acq_b);

	virtual unsigned int number() = 0;
	virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
	virtual void append_acquisition(ISMRMRD::Acquisition& acq) = 0;
	virtual void copy_parameters(const AcquisitionsContainer& ac) = 0;
	//virtual void write_parameters() = 0;
	virtual 
		boost::shared_ptr<AcquisitionsContainer> new_acquisitions_container() = 0;
	virtual int set_acquisition_data
		(int na, int nc, int ns, const float* re, const float* im) = 0;

	virtual void axpby(
		complex_float_t a, const aDataContainer<complex_float_t>& a_x,
		complex_float_t b, const aDataContainer<complex_float_t>& a_y);
	virtual complex_float_t dot(aDataContainer<complex_float_t>& dc);
	virtual float norm();
	float diff(AcquisitionsContainer& other);

	std::string parameters() const
	{
		return par_;
	}
	void set_parameters(std::string par)
	{
		par_ = par;
	}
	void set_ordered(bool ordered)
	{
		ordered_ = ordered;
	}
	bool undersampled() const
	{
		ISMRMRD::IsmrmrdHeader header;
		ISMRMRD::deserialize(par_.c_str(), header);
		ISMRMRD::Encoding e = header.encoding[0];
		return e.parallelImaging.is_present() &&
			e.parallelImaging().accelerationFactor.kspace_encoding_step_1 > 1;
	}
	int get_acquisitions_dimensions(size_t ptr_dim);

	void get_acquisitions_flags(unsigned int n, int* flags);

	unsigned int get_acquisitions_data(unsigned int slice, float* re, float* im);

	void order();

	bool ordered() const
	{
		return ordered_;
	}

	int* index()
	{
		return index_;
	}

	const int* index() const
	{
		return index_;
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
	void take_over(AcquisitionsContainer& ac)
	{
		AcquisitionsFile& af = (AcquisitionsFile&)ac;
		par_ = ac.parameters();
		if (index_)
			delete[] index_;
		int* index = ac.index();
		ordered_ = ac.ordered();
		if (ordered_ && index) {
			unsigned int n = number();
			index_ = new int[n];
			memcpy(index_, index, n*sizeof(int));
		}
		else
			index_ = 0;
		dataset_ = af.dataset_;
		if (own_file_) {
			Mutex mtx;
			mtx.lock();
			std::remove(filename_.c_str());
			mtx.unlock();
		}
		filename_ = af.filename_;
		own_file_ = af.own_file_;
		af.own_file_ = false;
	}
	virtual int set_acquisition_data
		(int na, int nc, int ns, const float* re, const float* im);
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
	virtual boost::shared_ptr<aDataContainer<complex_float_t> > new_data_container()
	{
		AcquisitionsFile* ptr_ac = acqs_scratch_file_(filename_);
		ptr_ac->set_parameters(par_);
		ptr_ac->write_parameters();
		boost::shared_ptr<aDataContainer<complex_float_t> > sptr_ac(ptr_ac);
		return sptr_ac;
	}
	virtual boost::shared_ptr<AcquisitionsContainer> new_acquisitions_container()
	{
		//boost::shared_ptr<AcquisitionsContainer> 
		//	sptr_ac(acqs_scratch_file_(filename_));
		//sptr_ac->set_parameters(par_);
		//sptr_ac->write_parameters();
		AcquisitionsFile* ptr_ac = acqs_scratch_file_(filename_);
		ptr_ac->set_parameters(par_);
		ptr_ac->write_parameters();
		boost::shared_ptr<AcquisitionsContainer> sptr_ac(ptr_ac);
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

class AcquisitionsVector : public AcquisitionsContainer {
public:
	virtual unsigned int number()
	{
		return (unsigned int)acqs_.size();
	}
	virtual unsigned int items()
	{
		return (unsigned int)acqs_.size();
	}
	virtual void append_acquisition(ISMRMRD::Acquisition& acq)
	{
		acqs_.push_back(boost::shared_ptr<ISMRMRD::Acquisition>
			(new ISMRMRD::Acquisition(acq)));
	}
	virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
	{
		int ind = index(num);
		acq = *acqs_[ind];
	}
	virtual void copy_parameters(const AcquisitionsContainer& ac)
	{
		par_ = ac.parameters();
	}
	virtual int set_acquisition_data
	(int na, int nc, int ns, const float* re, const float* im);
	virtual boost::shared_ptr<aDataContainer<complex_float_t> > new_data_container()
	{
		AcquisitionsVector* ptr_ac = new AcquisitionsVector();
		ptr_ac->set_parameters(par_);
		boost::shared_ptr<aDataContainer<complex_float_t> > sptr_ac(ptr_ac);
		return sptr_ac;
	}
	virtual boost::shared_ptr<AcquisitionsContainer> new_acquisitions_container()
	{
		AcquisitionsVector* ptr_ac = new AcquisitionsVector();
		ptr_ac->set_parameters(par_);
		boost::shared_ptr<AcquisitionsContainer> sptr_ac(ptr_ac);
		return sptr_ac;
	}
private:
	std::vector<boost::shared_ptr<ISMRMRD::Acquisition> > acqs_;
};

class ImagesContainer : public aDataContainer<complex_float_t> {
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
	virtual void get_images_data_as_float_array(float* data) const = 0;
	virtual void get_images_data_as_complex_array
		(float* re, float* im) const = 0;
	virtual void set_complex_images_data(const float* re, const float* im) = 0;
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
		complex_float_t a, const aDataContainer<complex_float_t>& a_x,
		complex_float_t b, const aDataContainer<complex_float_t>& a_y);
	virtual complex_float_t dot(aDataContainer<complex_float_t>& dc);
	virtual float norm();

	void get_image_data_as_cmplx_array
		(unsigned int im_num, float* re, float* im)
	{
		ImageWrap& iw = image_wrap(im_num);
		iw.get_cmplx_data(re, im);
	}

};

class ImagesVector : public ImagesContainer {
public:
	ImagesVector() : images_(), nimages_(0)
	{
	}
	ImagesVector(const ImagesVector& list, const char* attr, const char* target);
	ImagesVector(const ImagesVector& list, unsigned int inc = 1, unsigned int off = 0);
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
		return images_[im_num];
	}
	virtual boost::shared_ptr<const ImageWrap> sptr_image_wrap
	(unsigned int im_num) const
	{
		return images_[im_num];
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
	virtual void write(std::string filename, std::string groupname);
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
	virtual void get_images_data_as_float_array(float* data) const;
	virtual void get_images_data_as_complex_array(float* re, float* im) const;
	virtual void set_complex_images_data(const float* re, const float* im);
	virtual boost::shared_ptr<aDataContainer<complex_float_t> > new_data_container()
	{
		return boost::shared_ptr<aDataContainer<complex_float_t> >
			((aDataContainer<complex_float_t>*)new ImagesVector());
	}
	virtual boost::shared_ptr<ImagesContainer> new_images_container()
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesVector());
	}
	virtual boost::shared_ptr<ImagesContainer>
		clone(const char* attr, const char* target)
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesVector(*this, attr, target));
	}
	virtual boost::shared_ptr<ImagesContainer>
		clone(unsigned int inc = 1, unsigned int off = 0)
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesVector(*this, inc, off));
	}

private:
	std::vector<boost::shared_ptr<ImageWrap> > images_;
	int nimages_;
};

class CoilData {
public:
	virtual ~CoilData() {}
	virtual void get_dim(int* dim) const = 0;
	virtual void get_data(float* re, float* im) const = 0;
	virtual void set_data(const float* re, const float* im) = 0;
	virtual void get_data(complex_float_t* data) const = 0;
	virtual void set_data(const complex_float_t* data) = 0;
	virtual void get_data_abs(float* v) const = 0;
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
	virtual void get_data(float* re, float* im) const;
	virtual void set_data(const float* re, const float* im);
	virtual void get_data(complex_float_t* data) const
	{
		memcpy(data, img_.getDataPtr(), img_.getDataSize());
	}
	virtual void set_data(const complex_float_t* data)
	{
		memcpy(img_.getDataPtr(), data, img_.getDataSize());
	}
	virtual void get_data_abs(float* v) const;
private:
	ISMRMRD::Image < complex_float_t > img_;
};

class CoilDataContainer : public aDataContainer<complex_float_t> {
public:
	virtual float norm()
	{
		return 0.0;
	}
	virtual complex_float_t dot(aDataContainer<complex_float_t>& dc)
	{
		return complex_float_t(0.0, 0.0);
	}
	virtual void axpby(
		complex_float_t a, const aDataContainer<complex_float_t>& a_x,
		complex_float_t b, const aDataContainer<complex_float_t>& a_y)
	{
		return;
	}
	void get_dim(int slice, int* dim) //const
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.get_dim(dim);
	}
	void get_data(int slice, float* re, float* im) //const
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.get_data(re, im);
	}
	void set_data(int slice, float* re, float* im)
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
	void get_data_abs(int slice, float* v) //const
	{
		CoilData& ci = (CoilData&)(*this)(slice);
		ci.get_data_abs(v);
	}
	virtual void append(boost::shared_ptr<CoilData> sptr_csm) = 0;
	virtual CoilData& operator()(int slice) = 0;
	//virtual const CoilData& operator()(int slice) const = 0;
};

class CoilDataVector {
public:
	unsigned int items()
	{
		return (unsigned int)coil_data_.size();
	}
	CoilData& data(int slice) {
		return *coil_data_[slice];
	}
	virtual void append(boost::shared_ptr<CoilData> sptr_cd)
	{
		coil_data_.push_back(sptr_cd);
	}
private:
	std::vector< boost::shared_ptr<CoilData> > coil_data_;
};

class CoilImagesContainer : public CoilDataContainer {
public:
	virtual CoilData& operator()(int slice) = 0;
	virtual void compute(AcquisitionsContainer& ac);
	ISMRMRD::Encoding encoding() const
	{
		return encoding_;
	}
protected:
	ISMRMRD::Encoding encoding_;
};

class CoilImagesVector : public CoilImagesContainer, public CoilDataVector {
public:
	virtual boost::shared_ptr<aDataContainer<complex_float_t> > new_data_container()
	{
		return boost::shared_ptr<aDataContainer<complex_float_t> >
			((aDataContainer<complex_float_t>*)new CoilImagesVector());
	}
	virtual unsigned int items()
	{
		return CoilDataVector::items();
	}
	virtual CoilData& operator()(int slice)
	{
		return data(slice);
	}
	virtual void append(boost::shared_ptr<CoilData> sptr_cd)
	{
		CoilDataVector::append(sptr_cd);
	}
};

//template<class CoilDataType>
//class CoilSensitivitiesContainerTemplate : public CoilDataContainer {
class CoilSensitivitiesContainer : public CoilDataContainer {
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
		CoilImagesVector cis;
		cis.compute(ac);
		compute(cis);
	}

	virtual void compute(CoilImagesContainer& cis);

	void append_csm
		(int nx, int ny, int nz, int nc, const float* re, const float* im)
	{
		//CoilData* ptr_img = new CoilDataType(nx, ny, nz, nc);
		CoilData* ptr_img = new CoilDataAsCFImage(nx, ny, nz, nc);
		boost::shared_ptr<CoilData> sptr_img(ptr_img);
		ptr_img->set_data(re, im);
		append(sptr_img);
	}

protected:
	int csm_smoothness_;

private:
	void compute_csm_(
		ISMRMRD::NDArray<complex_float_t>& cm,
		ISMRMRD::NDArray<float>& img,
		ISMRMRD::NDArray<complex_float_t>& csm
	);

	float max_(int nx, int ny, float* u);
	void mask_noise_
	(int nx, int ny, float* u, float noise, int* mask);
	void cleanup_mask_(int nx, int ny, int* mask, int bg, int minsz, int ex);
	void smoothen_
	(int nx, int ny, int nz,
		complex_float_t* u, complex_float_t* v,
		int* obj_mask);
};

//typedef CoilSensitivitiesContainerTemplate<CoilDataAsCFImage> CoilSensitivitiesContainer;

class CoilSensitivitiesAsImages : public CoilSensitivitiesContainer, 
	public CoilDataVector {
public:
	CoilSensitivitiesAsImages()
	{
		csm_smoothness_ = 0;
	}
	CoilSensitivitiesAsImages(const char* file);

	virtual boost::shared_ptr<aDataContainer<complex_float_t> > new_data_container()
	{
		return boost::shared_ptr<aDataContainer<complex_float_t> >
			((aDataContainer<complex_float_t>*)new CoilSensitivitiesAsImages());
	}

	virtual unsigned int items()
	{
		return CoilDataVector::items();
	}
	virtual CoilData& operator()(int slice)
	{
		return data(slice);
	}
	virtual void append(boost::shared_ptr<CoilData> sptr_cd)
	{
		CoilDataVector::append(sptr_cd);
	}

};

#endif
