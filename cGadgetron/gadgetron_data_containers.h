#ifndef GADGETRON_DATA_CONTAINERS
#define GADGETRON_DATA_CONTAINERS

#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>

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
		t = z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, short& t)
	{
		t = z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, unsigned int& t)
	{
		t = z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, int& t)
	{
		t = z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, float& t)
	{
		t = z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, complex_float_t& t)
	{
		t = z;
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, double& t)
	{
		t = z.real();
	}
	template<typename T>
	static void convert_complex(std::complex<T> z, complex_double_t& t)
	{
		t = z;
	}
	static unsigned short abs(unsigned short v)
	{
		return v;
	}
	static short abs(short v)
	{
		return v > 0 ? v : -v;
	}
	static unsigned int abs(unsigned int v)
	{
		return v;
	}
	static int abs(int v)
	{
		return v > 0 ? v : -v;
	}
	static float abs(float v)
	{
		return v > 0 ? v : -v;
	}
	static double abs(double v)
	{
		return v > 0 ? v : -v;
	}
	static float abs(complex_float_t v)
	{
		return std::abs(v);
	}
	static double abs(complex_double_t v)
	{
		return std::abs(v);
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
	void get_dim(int* dim) const
	{
		IMAGE_PROCESSING_SWITCH_CONST(type_, get_dim_, ptr_, dim);
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

private:
	int type_;
	void* ptr_;

	ImageWrap& operator=(const ImageWrap& iw)
	{
		return *this;
	}

	template<typename T>
	void copy_(const ISMRMRD::Image<T>* ptr_im)
	{
		ptr_ = (void*)new ISMRMRD::Image<T>(*ptr_im);
	}

	template<typename T>
	void get_size_(const ISMRMRD::Image<T>* ptr_im, size_t& size) const
	{
		size = ptr_im->getDataSize();
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
		long long int n = im.getMatrixSizeX();
		n *= im.getMatrixSizeY();
		n *= im.getMatrixSizeZ();
		n *= im.getNumberOfChannels();
		const T* ptr = im.getDataPtr();
		for (long long int i = 0; i < n; i++)
			data[i] = xGadgetronUtilities::abs(ptr[i]);
	}

	template<typename T>
	void axpby_
		(const ISMRMRD::Image<T>* ptr_x, complex_double_t a, complex_double_t b)
	{
		ISMRMRD::Image<T>* ptr_y = (ISMRMRD::Image<T>*)ptr_;
		const T* i;
		T* j;
		long long ii = 0;
		long long int n = ptr_x->getMatrixSizeX();
		n *= ptr_x->getMatrixSizeY();
		n *= ptr_x->getMatrixSizeZ();
		n *= ptr_x->getNumberOfChannels();
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
		long long ii = 0;
		long long int n = ptr_im->getMatrixSizeX();
		n *= ptr_im->getMatrixSizeY();
		n *= ptr_im->getMatrixSizeZ();
		n *= ptr_im->getNumberOfChannels();
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
		long long ii = 0;
		long long int n = ptr->getMatrixSizeX();
		n *= ptr->getMatrixSizeY();
		n *= ptr->getMatrixSizeZ();
		n *= ptr->getNumberOfChannels();
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
		long long ii = 0;
		long long int n = ptr_im->getMatrixSizeX();
		n *= ptr_im->getMatrixSizeY();
		n *= ptr_im->getMatrixSizeZ();
		n *= ptr_im->getNumberOfChannels();
		*s = 0;
		for (i = ptr->getDataPtr(), j = ptr_im->getDataPtr(); ii < n; 
			i++, j++, ii++) {
			complex_float_t a = (complex_float_t)*i;
			complex_float_t b = (complex_float_t)*j;
			*s += std::abs(b - a);
		}
	}
};

class aDataContainer {
public:
	virtual ~aDataContainer() {}
	virtual boost::shared_ptr<aDataContainer> new_data_container() = 0;
	virtual double norm() = 0;
	virtual complex_double_t dot(aDataContainer& dc) = 0;
	virtual void axpby(
		complex_double_t a, const aDataContainer& a_x,
		complex_double_t b, const aDataContainer& a_y) = 0;
};

class AcquisitionsContainer : public aDataContainer {
public:
	std::string parameters() const
	{
		return par_;
	}
	void set_parameters(std::string par)
	{
		par_ = par;
	}
	boost::shared_ptr<ISMRMRD::NDArray<complex_float_t> > coils() const 
	{
		return coils_;
	}
	void set_coils(boost::shared_ptr<ISMRMRD::NDArray<complex_float_t> > coils)
	{
		coils_ = coils;
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

	virtual int number() = 0;
	virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
	virtual void append_acquisition(ISMRMRD::Acquisition& acq) = 0;
	virtual void copy_data(const AcquisitionsContainer& ac) = 0;
	virtual void write_data() = 0;
	virtual boost::shared_ptr<AcquisitionsContainer> new_acquisitions_container() = 0;

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
		for (int i = 0; i < n && i < m; i++) {
			y.get_acquisition(i, ay);
			x.get_acquisition(i, ax);
			AcquisitionsContainer::axpby(a, ax, b, ay);
			append_acquisition(ay);
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
		for (int i = 0; i < n && i < m; i++) {
			get_acquisition(i, a);
			other.get_acquisition(i, b);
			z += AcquisitionsContainer::dot(a, b);
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

protected:
	std::string par_;
	boost::shared_ptr<ISMRMRD::NDArray<complex_float_t> > coils_;
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
			coils_ = boost::shared_ptr<ISMRMRD::NDArray<complex_float_t> >
				(new ISMRMRD::NDArray<complex_float_t>);
			try {
				dataset_->readNDArray("csm", 0, *coils_);
				//ndim = coils_->getNDim();
			}
			catch (...) {
				std::cout << "WARNING: coil sensitivity map not found" << std::endl;
			}
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
	virtual int number()
	{
		Mutex mtx;
		mtx.lock();
		int na = dataset_->getNumberOfAcquisitions();
		mtx.unlock();
		return na;
	}
	virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
	{
		Mutex mtx;
		mtx.lock();
		dataset_->readAcquisition(num, acq);
		mtx.unlock();
	}
	virtual void append_acquisition(ISMRMRD::Acquisition& acq)
	{
		Mutex mtx;
		mtx.lock();
		dataset_->appendAcquisition(acq);
		mtx.unlock();
	}
	virtual void copy_data(const AcquisitionsContainer& ac) 
	{
		par_ = ac.parameters();
		coils_ = ac.coils();
		Mutex mtx;
		mtx.lock();
		dataset_->writeHeader(par_);
		dataset_->appendNDArray("csm", *coils_);
		mtx.unlock();
	}
	virtual void write_data() 
	{
		Mutex mtx;
		mtx.lock();
		dataset_->writeHeader(par_);
		dataset_->appendNDArray("csm", *coils_);
		mtx.unlock();
	}
	virtual boost::shared_ptr<aDataContainer> new_data_container()
	{
		AcquisitionsFile* ptr_ac = acqs_scratch_file_(filename_);
		ptr_ac->set_parameters(par_);
		ptr_ac->set_coils(coils_);
		ptr_ac->write_data();
		boost::shared_ptr<aDataContainer> sptr_ac(ptr_ac);
		return sptr_ac;
	}
	virtual boost::shared_ptr<AcquisitionsContainer> new_acquisitions_container()
	{
		boost::shared_ptr<AcquisitionsContainer> 
			sptr_ac(acqs_scratch_file_(filename_));
		sptr_ac->set_parameters(par_);
		sptr_ac->set_coils(coils_);
		sptr_ac->write_data();
		return sptr_ac;
	}

	void getPhantomAsFloat(ImageWrap& iw)
	{
		ISMRMRD::NDArray<complex_float_t> arr;
		dataset_->readNDArray("phantom", 0, arr);
		ISMRMRD::Image<float>* ptr_im =
			(ISMRMRD::Image<float>*)iw.ptr_image();
		complex_float_t* ia;
		float* ii;
		for (ia = arr.begin(), ii = ptr_im->begin(); ia != arr.end(); ia++, ii++)
			*ii = std::abs(*ia);
	}

	void getPhantomAsComplexFloat(ImageWrap& iw)
	{
		ISMRMRD::NDArray<complex_float_t> arr;
		dataset_->readNDArray("phantom", 0, arr);
		ISMRMRD::Image<complex_float_t>* ptr_im =
			(ISMRMRD::Image<complex_float_t>*)iw.ptr_image();
		complex_float_t* ia;
		complex_float_t* ii;
		for (ia = arr.begin(), ii = ptr_im->begin(); ia != arr.end(); ia++, ii++)
			*ii = std::abs(*ia);
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
	virtual int number()
	{
		return (int)acqs_.size();
	}
private:
	std::list<boost::shared_ptr<ISMRMRD::Acquisition> > acqs_;
};

class ImagesContainer : public aDataContainer {
public:
	virtual int number() const = 0;
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
	virtual void write(std::string filename, std::string groupname) = 0;
	virtual boost::shared_ptr<ImagesContainer> new_images_container() = 0;
	virtual boost::shared_ptr<ImagesContainer> clone() = 0;

	virtual void axpby(
		complex_double_t a, const aDataContainer& a_x,
		complex_double_t b, const aDataContainer& a_y)
	{
		ImagesContainer& x = (ImagesContainer&)a_x;
		ImagesContainer& y = (ImagesContainer&)a_y;
		ImageWrap w(x.image_wrap(0));
		complex_double_t zero(0.0, 0.0);
		complex_double_t one(1.0, 0.0);
		for (int i = 0; i < x.number() && i < y.number(); i++) {
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
		for (int i = 0; i < number() && i < ic.number(); i++) {
			const ImageWrap& u = image_wrap(i);
			const ImageWrap& v = ic.image_wrap(i);
			z += u.dot(v);
		}
		return z;
	}
	virtual double norm()
	{
		double r = 0;
		for (int i = 0; i < number(); i++) {
			const ImageWrap& u = image_wrap(i);
			double s = u.norm();
			r += s*s;
		}
		r = sqrt(r);
		return r;
	}
};

class ImagesList : public ImagesContainer {
public:
	ImagesList() : images_()
	{
	}
	ImagesList(const ImagesList& list)
	{
#ifdef MSVC
		std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::const_iterator i;
#endif
		for (i = list.images_.begin(); i != list.images_.end(); i++) {
			const boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			append(*sptr_iw);
		}
	}
	virtual int number() const 
	{
		return (int)images_.size();
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
#ifdef MSVC
		std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#endif
		unsigned int count = 0;
		for (i = images_.begin(); i != images_.end() && count < im_num; i++)
			count++;
		return *i;
	}
	virtual boost::shared_ptr<const ImageWrap> sptr_image_wrap
		(unsigned int im_num) const
	{
#ifdef MSVC
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
#ifdef MSVC
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
		if (im_num < 0 || im_num >= images_.size())
			dim[0] = dim[1] = dim[2] = 0;
		ImageWrap& iw = image_wrap(im_num);
		iw.get_dim(dim);
	}
	virtual void get_image_data_as_double_array(unsigned int im_num, double* data)
	{
		ImageWrap& iw = image_wrap(im_num);
		iw.get_data(data);
	}
	virtual boost::shared_ptr<aDataContainer> new_data_container()
	{
		return boost::shared_ptr<aDataContainer>((aDataContainer*)new ImagesList());
	}
	virtual boost::shared_ptr<ImagesContainer> new_images_container()
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesList());
	}
	virtual boost::shared_ptr<ImagesContainer> clone()
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesList(*this));
	}

private:
	std::list<boost::shared_ptr<ImageWrap> > images_;
};

#endif
