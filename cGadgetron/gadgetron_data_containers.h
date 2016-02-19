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
	virtual int number() = 0;
	virtual void getAcquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
	virtual void appendAcquisition(ISMRMRD::Acquisition& acq) = 0;
	virtual void writeHeader(const std::string& xml) = 0;

protected:
	std::string par_;
};

class AcquisitionsFile : public AcquisitionsContainer {
public:
	AcquisitionsFile(std::string filename, bool create = false, bool own = false)
	{
		own_file_ = own;
		filename_ = filename;
		Mutex mutex;
		boost::mutex& mtx = mutex();
		mtx.lock();
		dataset_ = boost::shared_ptr<ISMRMRD::Dataset>
			(new ISMRMRD::Dataset(filename.c_str(), "/dataset", create));
		if (!create)
			dataset_->readHeader(par_);
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
	virtual void writeHeader(const std::string& xml)
	{
		Mutex mutex;
		boost::mutex& mtx = mutex();
		mtx.lock();
		dataset_->writeHeader(xml);
		mtx.unlock();
		setParameters(xml);
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
		//if (type_ == ISMRMRD::ISMRMRD_USHORT)
		//	delete ((ISMRMRD::Image<unsigned short>*) ptr_);
		//else if (type_ == ISMRMRD::ISMRMRD_SHORT)
		//	delete (ISMRMRD::Image<short>*) ptr_;
		//else if (type_ == ISMRMRD::ISMRMRD_UINT)
		//	delete (ISMRMRD::Image<unsigned int>*) ptr_;
		//else if (type_ == ISMRMRD::ISMRMRD_INT)
		//	delete (ISMRMRD::Image<int>*) ptr_;
		//else if (type_ == ISMRMRD::ISMRMRD_FLOAT)
		//	delete (ISMRMRD::Image<float>*) ptr_;
		//else if (type_ == ISMRMRD::ISMRMRD_DOUBLE)
		//	delete (ISMRMRD::Image<double>*) ptr_;
		//else if (type_ == ISMRMRD::ISMRMRD_CXFLOAT)
		//	delete (ISMRMRD::Image< std::complex<float> >*) ptr_;
		//else if (type_ == ISMRMRD::ISMRMRD_CXDOUBLE)
		//	delete (ISMRMRD::Image< std::complex<double> >*) ptr_;
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
			//if (type == ISMRMRD::ISMRMRD_USHORT)
			//	writeImage(*(ISMRMRD::Image<unsigned short>*)ptr, dataset, mtx);
			//else if (type == ISMRMRD::ISMRMRD_SHORT)
			//	writeImage(*(ISMRMRD::Image<short>*)ptr, dataset, mtx);
			//else if (type == ISMRMRD::ISMRMRD_UINT)
			//	writeImage(*(ISMRMRD::Image<unsigned int>*)ptr, dataset, mtx);
			//else if (type == ISMRMRD::ISMRMRD_INT)
			//	writeImage(*(ISMRMRD::Image<int>*)ptr, dataset, mtx);
			//else if (type == ISMRMRD::ISMRMRD_FLOAT)
			//	writeImage(*(ISMRMRD::Image<float>*)ptr, dataset, mtx);
			//else if (type == ISMRMRD::ISMRMRD_DOUBLE)
			//	writeImage(*(ISMRMRD::Image<double>*)ptr, dataset, mtx);
			//else if (type == ISMRMRD::ISMRMRD_CXFLOAT)
			//	writeImage(*(ISMRMRD::Image< std::complex<float> >*)ptr, dataset, mtx);
			//else if (type == ISMRMRD::ISMRMRD_CXDOUBLE)
			//	writeImage(*(ISMRMRD::Image< std::complex<double> >*)ptr, dataset, mtx);
		}
	}
	virtual void getImageDimensions(unsigned int im_num, int* dim)
	{
		if (im_num < 0 || im_num >= images_.size())
			dim[0] = dim[1] = dim[2] = 0;
		ImageWrap& iw = imageWrap(im_num);
		int type = iw.type();
		void* ptr = iw.ptr_image();
		if (type == ISMRMRD::ISMRMRD_USHORT)
			getImageDim(*(ISMRMRD::Image<unsigned short>*)ptr, dim);
		else if (type == ISMRMRD::ISMRMRD_SHORT)
			getImageDim(*(ISMRMRD::Image<short>*)ptr, dim);
		else if (type == ISMRMRD::ISMRMRD_UINT)
			getImageDim(*(ISMRMRD::Image<unsigned int>*)ptr, dim);
		else if (type == ISMRMRD::ISMRMRD_INT)
			getImageDim(*(ISMRMRD::Image<int>*)ptr, dim);
		else if (type == ISMRMRD::ISMRMRD_FLOAT)
			getImageDim(*(ISMRMRD::Image<float>*)ptr, dim);
		else if (type == ISMRMRD::ISMRMRD_DOUBLE)
			getImageDim(*(ISMRMRD::Image<double>*)ptr, dim);
		else if (type == ISMRMRD::ISMRMRD_CXFLOAT)
			getImageDim(*(ISMRMRD::Image< std::complex<float> >*)ptr, dim);
		else if (type == ISMRMRD::ISMRMRD_CXDOUBLE)
			getImageDim(*(ISMRMRD::Image< std::complex<double> >*)ptr, dim);
	}
	virtual void getImageDataAsDoubleArray(unsigned int im_num, double* data)
	{
		ImageWrap& iw = imageWrap(im_num);
		int type = iw.type();
		void* ptr = iw.ptr_image();
		if (type == ISMRMRD::ISMRMRD_USHORT)
			getImageUnsignedData(*(ISMRMRD::Image<unsigned short>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_SHORT)
			getImageData(*(ISMRMRD::Image<short>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_UINT)
			getImageUnsignedData(*(ISMRMRD::Image<unsigned int>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_INT)
			getImageData(*(ISMRMRD::Image<int>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_FLOAT)
			getImageData(*(ISMRMRD::Image<float>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_DOUBLE)
			getImageData(*(ISMRMRD::Image<double>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_CXFLOAT)
			getImageComplexData(*(ISMRMRD::Image< std::complex<float> >*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_CXDOUBLE)
			getImageComplexData(*(ISMRMRD::Image< std::complex<double> >*)ptr, data);
	}

private:

	std::list<boost::shared_ptr<ImageWrap> > images_;

	template<typename T>
	void writeImage
		(ISMRMRD::Image<T>* ptr_im, ISMRMRD::Dataset& dataset, boost::mutex& mtx)
		//(ISMRMRD::Image<T>& im, ISMRMRD::Dataset& dataset, boost::mutex& mtx)
	{
		ISMRMRD::Image<T>& im = *ptr_im;
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
	void getImageDim(ISMRMRD::Image<T>& im, int* dim) 
	{
		dim[0] = im.getMatrixSizeX();
		dim[1] = im.getMatrixSizeY();
		dim[2] = im.getMatrixSizeZ();
	}

#define ABS(X) ((X < 0) ? -X : X)

	template<typename T>
	void getImageData(ISMRMRD::Image<T>& im, double* data) 
	{
		long long int n = im.getMatrixSizeX();
		n *= im.getMatrixSizeY();
		n *= im.getMatrixSizeZ();
		T* ptr = im.getDataPtr();
		for (long long int i = 0; i < n; i++)
			data[i] = ABS(ptr[i]);
	}
	template<typename T>
	void getImageUnsignedData(ISMRMRD::Image<T>& im, double* data) 
	{
		long long int n = im.getMatrixSizeX();
		n *= im.getMatrixSizeY();
		n *= im.getMatrixSizeZ();
		T* ptr = im.getDataPtr();
		for (long long int i = 0; i < n; i++)
			data[i] = ptr[i];
	}
	template<typename T>
	void getImageComplexData(ISMRMRD::Image< std::complex<T> >& im, double* data)
	{
		long long int n = im.getMatrixSizeX();
		n *= im.getMatrixSizeY();
		n *= im.getMatrixSizeZ();
		std::complex<T>* ptr = im.getDataPtr();
		for (long long int i = 0; i < n; i++)
			data[i] = std::abs(ptr[i]);
	}
};

#endif
