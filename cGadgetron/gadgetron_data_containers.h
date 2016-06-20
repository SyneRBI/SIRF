#ifndef GADGETRON_DATA_CONTAINERS
#define GADGETRON_DATA_CONTAINERS

#include <complex>

#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>
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
	ISMRMRD::ImageHeader* ptr_head()
	{
		ISMRMRD::ImageHeader** h;
		IMAGE_PROCESSING_SWITCH(type_, get_head_ptr_, ptr_, h);
		return *h;
	}
	void set_imtype(ISMRMRD::ISMRMRD_ImageTypes imtype)
	{
		IMAGE_PROCESSING_SWITCH(type_, set_imtype_, ptr_, imtype);
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
		const CFImage& img = *(const CFImage*)ptr_;
		long long int n = img.getMatrixSizeX();
		n *= img.getMatrixSizeY();
		n *= img.getMatrixSizeZ();
		n *= img.getNumberOfChannels();
		const complex_float_t* ptr = img.getDataPtr();
		for (long long int i = 0; i < n; i++) {
			complex_float_t z = ptr[i];
			re[i] = std::real(z);
			im[i] = std::imag(z);
		}
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
	virtual int items() = 0;
	virtual double norm() = 0;
	virtual complex_double_t dot(aDataContainer& dc) = 0;
	virtual void axpby(
		complex_double_t a, const aDataContainer& a_x,
		complex_double_t b, const aDataContainer& a_y) = 0;
};

class AcquisitionsContainer : public aDataContainer {
public:
	AcquisitionsContainer() : index_(0) {}
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

	void get_acquisitions_dimensions(size_t ptr_dim)
	{
		ISMRMRD::Acquisition acq;
		int* dim = (int*)ptr_dim;

		int y = 0;
		for (;;){
			get_acquisition(y, acq);
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
				break;
			y++;
		}
		int ny = 0;
		for (;;) {
			get_acquisition(y, acq);
			y++;
			ny++;
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
				break;
		}

		dim[0] = acq.number_of_samples();
		dim[1] = ny; // e.reconSpace.matrixSize.y;
		dim[2] = acq.active_channels();
	}
	void get_acquisitions_data(unsigned int slice, double* re, double* im)
	{
		ISMRMRD::Acquisition acq;
		int* dim = new int[3];
		size_t ptr_dim = (size_t)dim;
		get_acquisitions_dimensions(ptr_dim);
		unsigned int ny = dim[1]; //e.reconSpace.matrixSize.y;
		delete[] dim;
		int y = 0;
		for (;;){
			get_acquisition(y + ny*slice, acq);
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
				break;
			y++;
		}
		for (;;) {
			get_acquisition(y + ny*slice, acq);
			unsigned int nc = acq.active_channels();
			unsigned int ns = acq.number_of_samples();
			for (size_t c = 0; c < nc; c++) {
				for (size_t s = 0; s < ns; s++) {
					complex_float_t z = acq.data(s, c);
					re[s + ns*(y + ny*c)] = std::real(z);
					im[s + ns*(y + ny*c)] = std::imag(z);
				}
			}
			y++;
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
				break;
		}
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
		//for (int i = 0; i < n && i < m; i++) {
		for (int i = 0, j = 0; i < n && j < m;) {
			get_acquisition(i, a);
			if (TO_BE_IGNORED(a)) {
				std::cout << i << " ignored (a)\n";
				i++;
				continue;
			}
			other.get_acquisition(j, b);
			if (TO_BE_IGNORED(b)) {
				std::cout << j << " ignored (b)\n";
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
				std::cout << i << " ignored (norm)\n";
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
			//std::cout << i << ' ' << j << '\n';
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
			//std::cout << i << ' ' << j << '\n';
			for (int k = i; k < j; k++)
				ts.push_back(triple(t[k]));
			std::sort(ts.begin(), ts.end(),
				[](triple a, triple b) { return b.third > a.third; });
			for (int k = i; k < j; k++)
				index_[k] = ts[k - i].ind;
			i = j;
		}
	}
	bool ordered() const
	{
		return (bool)index_;
	}
	int index(int i)
	{
		if (index_ && i >= 0 && i < number())
			return index_[i];
		else
			return i;
	}

protected:
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
	virtual int items()
	{
		Mutex mtx;
		mtx.lock();
		int na = dataset_->getNumberOfAcquisitions();
		mtx.unlock();
		return na;
	}
	virtual int number()
	{
		return items();
	}
	virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
	{
		int ind = index(num);
		//if (ordered())
		//	std::cout << num << ' ' << ind << '\n';
		Mutex mtx;
		mtx.lock();
		dataset_->readAcquisition(ind, acq);
		//dataset_->readAcquisition(num, acq);
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
	virtual int number() = 0;
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
	virtual void write(std::string filename, std::string groupname) = 0;
	virtual boost::shared_ptr<ImagesContainer> new_images_container() = 0;
	virtual boost::shared_ptr<ImagesContainer> 
		clone(unsigned int inc = 1, unsigned int off = 0) = 0;

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

	void get_image_data_as_cmplx_array
		(unsigned int im_num, double* re, double* im)
	{
		ImageWrap& iw = image_wrap(im_num);
		iw.get_cmplx_data(re, im);
	}

	void set_image_to_real_conversion(int type)
	{
		for (int i = 0; i < number(); i++) {
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
	ImagesList(const ImagesList& list, unsigned int inc = 1, unsigned int off = 0)
	{
		int n = 0;
		int j = 0;
#ifdef MSVC
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
	virtual int items()
	{
		return (int)images_.size();
	}
	virtual int number()
	{
		return (int)images_.size();
	}
	virtual int types()
	{
		if (nimages_ > 0)
			return images_.size() / nimages_;
		else
			return 1;
	}
	virtual void count(int i)
	{
		if (i > nimages_ - 1)
			nimages_ = i + 1;
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
		for (i = images_.begin(); 
			i != images_.end() && count < im_num && count < images_.size() - 1; i++)
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
		if (im_num >= images_.size())
			dim[0] = dim[1] = dim[2] = 0;
		ImageWrap& iw = image_wrap(im_num);
		iw.get_dim(dim);
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
	virtual boost::shared_ptr<aDataContainer> new_data_container()
	{
		return boost::shared_ptr<aDataContainer>((aDataContainer*)new ImagesList());
	}
	virtual boost::shared_ptr<ImagesContainer> new_images_container()
	{
		return boost::shared_ptr<ImagesContainer>(new ImagesList());
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

class CoilSensitivityMap {
public:
	virtual ~CoilSensitivityMap() {}
	virtual complex_float_t& operator()(int x, int y, int z, int c) = 0;
};

class CSMAsCFImage : public CoilSensitivityMap {
public:
	CSMAsCFImage(uint16_t nx = 0, uint16_t ny = 1, uint16_t nz = 1, uint16_t nc = 1) :
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
	void get_dim(int* dim) const
	{
		dim[0] = img_.getMatrixSizeX();
		dim[1] = img_.getMatrixSizeY();
		dim[2] = img_.getMatrixSizeZ();
		dim[3] = img_.getNumberOfChannels();
	}
	void get_data(double* re, double* im) const
	{
		long long int n = img_.getMatrixSizeX();
		n *= img_.getMatrixSizeY();
		n *= img_.getMatrixSizeZ();
		n *= img_.getNumberOfChannels();
		const complex_float_t* ptr = img_.getDataPtr();
		for (long long int i = 0; i < n; i++) {
			complex_float_t z = ptr[i];
			re[i] = std::real(z);
			im[i] = std::imag(z);
		}
	}
	void get_data_abs(double* v) const
	{
		long long int n = img_.getMatrixSizeX();
		n *= img_.getMatrixSizeY();
		n *= img_.getMatrixSizeZ();
		n *= img_.getNumberOfChannels();
		const complex_float_t* ptr = img_.getDataPtr();
		for (long long int i = 0; i < n; i++) {
			complex_float_t z = ptr[i];
			v[i] = std::abs(z);
		}
	}
private:
	ISMRMRD::Image < complex_float_t > img_;
};

class CoilSensitivitiesContainer : public aDataContainer {
public:
	void set_csm_smoothness(int s)
	{
		csm_smoothness_ = s;
	}
	virtual void get_dim(int slice, int* dim) = 0;
	virtual void get_data(int slice, double* re, double* im) = 0;
	virtual void get_data_abs(int slice, double* v) = 0;
	virtual void append(boost::shared_ptr<CoilSensitivityMap> sptr_csm) = 0;
	virtual void compute(AcquisitionsContainer& ac) = 0;
	virtual CoilSensitivityMap& operator()(int slice) = 0;
protected:
	int csm_smoothness_;
	std::list< boost::shared_ptr<CoilSensitivityMap> > csms_;
};

class CoilSensitivitiesAsImages : public CoilSensitivitiesContainer {
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
			boost::shared_ptr<CoilSensitivityMap> sptr_img(new CSMAsCFImage);
			mtx.lock();
			CFImage& csm = (*(CSMAsCFImage*)sptr_img.get()).image();
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
	virtual int items()
	{
		return (int)csms_.size();
	}

	virtual void append(boost::shared_ptr<CoilSensitivityMap> sptr_csm)
	{
		csms_.push_back(sptr_csm);
	}
	virtual void get_dim(int slice, int* dim)
	{
		CSMAsCFImage& csm = (CSMAsCFImage&)(*this)(slice);
		csm.get_dim(dim);
	}
	virtual void get_data(int slice, double* re, double* im)
	{
		CSMAsCFImage& csm = (CSMAsCFImage&)(*this)(slice);
		csm.get_data(re, im);
	}
	virtual void get_data_abs(int slice, double* v)
	{
		CSMAsCFImage& csm = (CSMAsCFImage&)(*this)(slice);
		csm.get_data_abs(v);
	}
	virtual void compute(AcquisitionsContainer& ac)
	{
		//if (!ac.ordered())
		//	ac.order();

		std::string par;
		ISMRMRD::IsmrmrdHeader header;
		ISMRMRD::Acquisition acq;
		par = ac.parameters();
		ISMRMRD::deserialize(par.c_str(), header);
		ac.get_acquisition(0, acq);

		ISMRMRD::Encoding e = header.encoding[0];
		//unsigned int readout = e.encodedSpace.matrixSize.x;
		unsigned int nx = e.reconSpace.matrixSize.x;
		unsigned int ny = e.reconSpace.matrixSize.y;
		unsigned int nc = acq.active_channels();
		unsigned int readout = acq.number_of_samples();

		//std::cout << readout << ' ' << acq.number_of_samples() << '\n';

		std::vector<size_t> cm_dims;
		cm_dims.push_back(readout);
		cm_dims.push_back(ny);
		cm_dims.push_back(nc);

		ISMRMRD::NDArray<complex_float_t> cm(cm_dims);

		std::vector<size_t> img_dims;
		img_dims.push_back(nx);
		img_dims.push_back(ny);
		ISMRMRD::NDArray<float> img(img_dims);

		int nmap = 0;
		std::cout << "map ";

		for (int na = 0; na < ac.number();) {

			std::cout << ++nmap << ' ' << std::flush;

			memset(cm.getDataPtr(), 0, cm.getDataSize());
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
				if (!e.parallelImaging.is_present() || 
					acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) ||
					acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING)) {
					for (size_t c = 0; c < nc; c++) {
						for (size_t s = 0; s < readout; s++) {
							cm(s, yy, c) = acq.data(s, c);
						}
					}
				}
				y++;
				if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
					break;
			}
			na += y;
			
			ifft2c(cm);

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
			//for (long long int i = 0; i < nx*ny*nc; i++)
			//	mask[i] = 1;

			ISMRMRD::NDArray<complex_float_t> w(cm0);
			//for (int i = 0; i < csm_smoothness_; i++)
				//smoothen_(cm0, w);

			float* ptr_img = img.getDataPtr();
			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++) {
					double r = 0.0;
					for (unsigned int c = 0; c < nc; c++) {
						float s = std::abs(cm0(x, y, c));
						r += s*s;
					}
					img(x, y) = std::sqrt(r);
				}
			}

			int* edge_mask = new int[nx*ny];			
			memset(edge_mask, 0, nx*ny*sizeof(int));
			float* weight = new float[nx*ny];
			find_edges_(nx, ny, ptr_img, edge_mask, weight);
			cleanup_mask_(nx, ny, edge_mask, 0, 2, 4);
			cleanup_mask_(nx, ny, edge_mask, 0, 4, 2);
			cleanup_mask_(nx, ny, edge_mask, 0, 4, 2);

			float noise = max_(5, 5, ptr_img);
			//std::cout << "noise: " << noise << '\n';
			mask_noise_(nx, ny, ptr_img, noise, object_mask);
			cleanup_mask_(nx, ny, object_mask, 0, 2, 0);
			cleanup_mask_(nx, ny, object_mask, 0, 3, 0);
			cleanup_mask_(nx, ny, object_mask, 0, 4, 0);
			for (int i = 0; i < csm_smoothness_; i++)
				smoothen_(nx, ny, nc, cm0.getDataPtr(), w.getDataPtr(), 
				object_mask, edge_mask);

			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++) {
					double r = 0.0;
					for (unsigned int c = 0; c < nc; c++) {
						float s = std::abs(cm0(x, y, c));
						r += s*s;
					}
					img(x, y) = std::sqrt(r);
				}
			}

			boost::shared_ptr<CoilSensitivityMap>
				sptr_img(new CSMAsCFImage(nx, ny, 1, nc));
			CFImage& csm = (*(CSMAsCFImage*)sptr_img.get()).image();

			for (unsigned int y = 0, i = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++, i++) {
					double r = img(x, y);
					float s;
					if (r != 0.0)
						s = 1.0 / r; //std::sqrt(r);
					else
						s = 0.0;
					complex_float_t z(s, 0.0);
					for (unsigned int c = 0; c < nc; c++) {
						//csm(x, y, 0, c) = weight[i];
						//csm(x, y, 0, c) = object_mask[i];
						//csm(x, y, 0, c) = edge_mask[i];
						//csm(x, y, 0, c) = img(x, y);
						csm(x, y, 0, c) = cm0(x, y, c) * z;
					}
				}
			}

			//float* u = new float[nx*ny];
			//for (unsigned int c = 0; c < nc; c++) {
			//	for (unsigned int y = 0; y < ny; y++) {
			//		for (unsigned int x = 0; x < nx; x++) {
			//			csm(x, y, 0, c) = img(x, y);
			//		}
			//	}
			//	smoothen_(nx, ny, ptr_img, u, object_mask, edge_mask);
			//}
			//delete[] u;

			//for (int i = 0; i < csm_smoothness_; i++)
			//	smoothen_(nx, ny, nc, csm.getDataPtr(), w.getDataPtr(), edge_mask);

			append(sptr_img);

			delete[] object_mask;
			delete[] edge_mask;
			delete[] weight;

		}
		std::cout << '\n';
	}

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
	virtual CoilSensitivityMap& operator()(int slice)
	{
#ifdef MSVC
		std::list<boost::shared_ptr<CoilSensitivityMap> >::const_iterator i;
#else
		typename std::list<boost::shared_ptr<CoilSensitivityMap> >::const_iterator i;
#endif
		unsigned int count = 0;
		for (i = csms_.begin(); 
			i != csms_.end() && count < slice && count < csms_.size() - 1; i++)
			count++;
		return **i;
	}

private:

	void smoothen_(ISMRMRD::NDArray<complex_float_t>& u, 
		ISMRMRD::NDArray<complex_float_t>& v)
	{
		const complex_float_t TWO(2.0, 0.0);
		const complex_float_t SIXTEEN(16.0, 0.0);
		const size_t *dims = u.getDims();
		int nx = dims[0];
		int ny = dims[1];
		int nz = dims[2];
		for (int iz = 0; iz < nz; iz++)
			for (int iy = 0; iy < ny; iy++)
				for (int ix = 0; ix < nx; ix++)
					if (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1)
						v(ix, iy, iz) = u(ix, iy, iz);
					else
						v(ix, iy, iz) = u(ix, iy, iz) / TWO + (
						u(ix - 1, iy - 1, iz) + u(ix, iy - 1, iz) + u(ix + 1, iy - 1, iz) +
						u(ix - 1, iy    , iz)                     + u(ix + 1, iy    , iz) +
						u(ix - 1, iy + 1, iz) + u(ix, iy + 1, iz) + u(ix + 1, iy + 1, iz)
						) / SIXTEEN;
		u = v;
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
		(int nx, int ny, int nz, complex_float_t* u, complex_float_t* v, 
		int* obj_mask, int* edge_mask)
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
	void smoothen_
		(int nx, int ny, float* u, float* v,
		int* obj_mask, int* edge_mask)
	{
		for (int iy = 0, i = 0; iy < ny; iy++)
			for (int ix = 0; ix < nx; ix++, i++) {
				if (edge_mask[i]) {
					v[i] = u[i];
					continue;
				}
				int n = 0;
				float s = 0.0;
				for (int jy = -1; jy <= 1; jy++)
					for (int jx = -1; jx <= 1; jx++) {
						if (ix + jx < 0 || ix + jx >= nx)
							continue;
						if (iy + jy < 0 || iy + jy >= ny)
							continue;
						int j = i + jx + jy*nx;
						if (i != j && !edge_mask[j]) {
							n++;
							s += u[j];
						}
					}
				//v[i] = (u[i] + s / 8.0) / 2;
				if (n > 0)
					v[i] = (u[i] + s / n) / 2;
				else
					v[i] = u[i];
			}
		memcpy(u, v, nx*ny*sizeof(float));
	}
	float weight_(float x)
	{
		if (x*x > 1)
			return 0;
		float y = 1 - x*x;
		return y*y;
	}
	void find_edges_(int nx, int ny, float* u, int* mask, float* weight)
	{
		int ng = (nx - 2)*(ny - 2);
		float* grad = new float[ng];
		float max_grad = 0.0;
		for (int iy = 1, i = 0; iy < ny - 1; iy++)
			for (int ix = 1; ix < nx - 1; ix++, i++) {
				int j = ix + iy*nx;
				float dx = fabs(u[j + 1] - u[j - 1]);
				float dy = fabs(u[j + nx] - u[j - nx]);
				//float dx = fabs(u[j + 1] - u[j]);
				//float dy = fabs(u[j + nx] - u[j]);
				float d = std::sqrt(dx*dx + dy*dy);
				if (d > max_grad)
					max_grad = d;
				grad[i] = d;
			}
		const float AH = 5;
		int nh = ng / AH;
		int* h = new int[nh];
		memset(h, 0, nh*sizeof(int));
		for (int i = 0; i < ng; i++) {
			int j = (nh - 1) * (grad[i] / max_grad);
			h[j]++;
		}
		std::ofstream out;
		//out.open("h.txt");
		//for (int i = 0; i < nh; i++)
		//	out << h[i] << '\n';
		//out.close();
		float cutoff = 0.0;
		float step = max_grad / nh;
		for (int i = 10; i < nh; i++, cutoff += step)
			if (h[i] < 1)
				break;
		//std::cout << "\nmax grad: " << max_grad << '\n';
		//std::cout << "cutoff: " << cutoff << '\n';
		for (int ix = 0; ix < nx; ix++) {
			weight[ix] = 1.0;
			weight[ix + nx*(ny - 1)] = 1.0;
		}
		for (int iy = 0; iy < ny; iy++) {
			weight[nx*iy] = 1.0;
			weight[nx*iy + nx - 1] = 1.0;
		}
		for (int iy = 1, i = 0; iy < ny - 1; iy++)
			for (int ix = 1; ix < nx - 1; ix++, i++) {
				int j = ix + iy*nx;
				weight[j] = weight_(grad[i] / cutoff);
				if (grad[i] > cutoff) {
					mask[j] = 1;
				}
			}

		delete[] grad;
		delete[] h;
	}
};

#endif
