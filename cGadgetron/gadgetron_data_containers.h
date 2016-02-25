#ifndef GADGETRON_DATA_CONTAINERS
#define GADGETRON_DATA_CONTAINERS

#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

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

class aDataContainer {
public:
	virtual ~aDataContainer() {}
};

class AcquisitionsContainer : public aDataContainer {
public:
	std::string parameters() const
	{
		return par_;
	}
	void setParameters(std::string par)
	{
		par_ = par;
	}
	boost::shared_ptr<ISMRMRD::NDArray<complex_float_t> > coils() const {
		return coils_;
	}
	virtual int number() = 0;
	virtual void getAcquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
	virtual void appendAcquisition(ISMRMRD::Acquisition& acq) = 0;
	virtual void copyData(const AcquisitionsContainer& ac) = 0;

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
		Mutex mutex;
		boost::mutex& mtx = mutex();
		mtx.lock();
		dataset_ = boost::shared_ptr<ISMRMRD::Dataset>
			(new ISMRMRD::Dataset(filename.c_str(), "/dataset", create_file));
		if (!create_file) {
			dataset_->readHeader(par_);
			coils_ = boost::shared_ptr<ISMRMRD::NDArray<complex_float_t> >
				(new ISMRMRD::NDArray<complex_float_t>);
			dataset_->readNDArray("csm", 0, *coils_);
		}
		mtx.unlock();
	}
	~AcquisitionsFile() {
		dataset_.reset();
		if (own_file_) {
			//std::cout << "removing " << filename_.c_str() << std::endl;
			std::remove(filename_.c_str());
		}

	}
	virtual int number()
	{
		Mutex mutex;
		boost::mutex& mtx = mutex();
		mtx.lock();
		int na = dataset_->getNumberOfAcquisitions();
		mtx.unlock();
		return na;
	}
	virtual void getAcquisition(unsigned int num, ISMRMRD::Acquisition& acq)
	{
		Mutex mutex;
		boost::mutex& mtx = mutex();
		mtx.lock();
		dataset_->readAcquisition(num, acq);
		mtx.unlock();
	}
	virtual void appendAcquisition(ISMRMRD::Acquisition& acq)
	{
		Mutex mutex;
		boost::mutex& mtx = mutex();
		mtx.lock();
		dataset_->appendAcquisition(acq);
		mtx.unlock();
	}
	virtual void copyData(const AcquisitionsContainer& ac) {
		par_ = ac.parameters();
		coils_ = ac.coils();
		Mutex mutex;
		boost::mutex& mtx = mutex();
		mtx.lock();
		dataset_->writeHeader(par_);
		dataset_->appendNDArray("csm", *coils_);
		mtx.unlock();
	}

private:
	bool own_file_;
	std::string filename_;
	boost::shared_ptr<ISMRMRD::Dataset> dataset_;
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

class ImageWrap {
public:
	ImageWrap(uint16_t type, void* ptr_im) 
	{
		type_ = type;
		ptr_ = ptr_im;
	}
	~ImageWrap() 
	{
		IMAGE_PROCESSING_SWITCH(type_, delete, ptr_);
	}
	int type() 
	{
		return type_;
	}
	void* ptr_image() 
	{
		return ptr_;
	}
private:
	int type_;
	void* ptr_;
};

class ImagesContainer : public aDataContainer {
public:
	virtual int number() const = 0;
	virtual ImageWrap& imageWrap(unsigned int im_num) = 0;
	virtual void append(int image_data_type, void* ptr_image) = 0;
	virtual void getImageDimensions(unsigned int im_num, int* dim) = 0;
	virtual void getImageDataAsDoubleArray(unsigned int im_num, double* data) = 0;
	virtual void write(std::string filename, std::string groupname) = 0;
};

class ImagesList : public ImagesContainer {
public:
	virtual int number() const 
	{
		return (int)images_.size();
	}
	virtual void append(int image_data_type, void* ptr_image) 
	{
		images_.push_back(boost::shared_ptr<ImageWrap>
			(new ImageWrap(image_data_type, ptr_image)));
	}
	virtual ImageWrap& imageWrap(unsigned int im_num)
	{
#ifdef MSVC
		std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#endif
		unsigned int count = 0;
		for (i = images_.begin(); i != images_.end() && count < im_num; i++)
			count++;
		boost::shared_ptr<ImageWrap>& sptr_iw = *i;
		return *sptr_iw;
	}
	virtual void write(std::string filename, std::string groupname)
	{
		if (images_.size() < 1)
			return;
		Mutex mutex;
		boost::mutex& mtx = mutex();
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
			int type = iw.type();
			void* ptr = iw.ptr_image();
			IMAGE_PROCESSING_SWITCH(type, writeImage, ptr, dataset, mtx);
		}
	}
	virtual void getImageDimensions(unsigned int im_num, int* dim)
	{
		if (im_num < 0 || im_num >= images_.size())
			dim[0] = dim[1] = dim[2] = 0;
		ImageWrap& iw = imageWrap(im_num);
		int type = iw.type();
		void* ptr = iw.ptr_image();
		IMAGE_PROCESSING_SWITCH(type, getImageDim, ptr, dim);
	}
	virtual void getImageDataAsDoubleArray(unsigned int im_num, double* data)
	{
		ImageWrap& iw = imageWrap(im_num);
		int type = iw.type();
		void* ptr = iw.ptr_image();
		IMAGE_PROCESSING_SWITCH(type, getImageData, ptr, data);
	}

private:

	std::list<boost::shared_ptr<ImageWrap> > images_;

	template<typename T>
	void writeImage
		(const ISMRMRD::Image<T>* ptr_im, ISMRMRD::Dataset& dataset, 
		boost::mutex& mtx)
	{
		const ISMRMRD::Image<T>& im = *ptr_im;
		std::stringstream ss;
		ss << "image_" << im.getHead().image_series_index;
		std::string image_varname = ss.str();
		{
			mtx.lock();
			dataset.appendImage(image_varname, im);
			mtx.unlock();
		}
	}
	template<typename T>
	void getImageDim(const ISMRMRD::Image<T>* ptr_im, int* dim)
	{
		const ISMRMRD::Image<T>& im = *ptr_im;
		dim[0] = im.getMatrixSizeX();
		dim[1] = im.getMatrixSizeY();
		dim[2] = im.getMatrixSizeZ();
	}

	unsigned short myabs(unsigned short v)
	{
		return v;
	}
	short myabs(short v)
	{
		return v > 0 ? v : -v;
	}
	unsigned int myabs(unsigned int v)
	{
		return v;
	}
	int myabs(int v)
	{
		return v > 0 ? v : -v;
	}
	float myabs(float v)
	{
		return v > 0 ? v : -v;
	}
	double myabs(double v)
	{
		return v > 0 ? v : -v;
	}
	float myabs(std::complex<float> v)
	{
		return std::abs(v);
	}
	double myabs(std::complex<double> v)
	{
		return std::abs(v);
	}

	template<typename T>
	void getImageData(const ISMRMRD::Image<T>* ptr_im, double* data)
	{
		const ISMRMRD::Image<T>& im = *ptr_im;
		long long int n = im.getMatrixSizeX();
		n *= im.getMatrixSizeY();
		n *= im.getMatrixSizeZ();
		const T* ptr = im.getDataPtr();
		for (long long int i = 0; i < n; i++)
			data[i] = myabs(ptr[i]);
	}
};

#endif
