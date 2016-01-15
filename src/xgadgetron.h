#ifndef GADGETRON_EXTENSIONS
#define GADGETRON_EXTENSIONS

#include <cmath>

#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>

#include "gadgetron_client.h"

class Mutex {
public:
	Mutex() {
		init_();
	}
	boost::mutex& operator()() {
		return *sptr_mutex_.get();
	}
	boost::shared_ptr<boost::mutex> sptr() {
		return sptr_mutex_;
	}
private:
	static boost::shared_ptr<boost::mutex> sptr_mutex_;
	static void init_() {
		static bool initialized = false;
		if (!initialized) {
			sptr_mutex_ = boost::shared_ptr<boost::mutex>(new boost::mutex);
			initialized = true;
		}
	}
};

class GTConnector {
public:
	GTConnector() {
		sptr_mtx_ = boost::shared_ptr<Mutex>(new Mutex);
		sptr_con_ = boost::shared_ptr<GadgetronClientConnector>
			(new GadgetronClientConnector);
	}
	GadgetronClientConnector& operator()() {
		return *sptr_con_.get();
	}
	boost::shared_ptr<GadgetronClientConnector> sptr() {
		return sptr_con_;
	}
	boost::mutex& mutex() {
		Mutex& mtx = *sptr_mtx_.get();
		return mtx();
	}
private:
	boost::shared_ptr<Mutex> sptr_mtx_;
	boost::shared_ptr<GadgetronClientConnector> sptr_con_;
};

class ImagesList {
public:
	std::list<boost::shared_ptr<ImageWrap> >& images() {
		return images_;
	}
	const std::list<boost::shared_ptr<ImageWrap> >& images() const {
		return images_;
	}
	int size() const {
		return images_.size();
	}
	void write
		(std::string filename, std::string groupname, GTConnector& conn)
	{
		if (images_.size() < 1)
			return;
		boost::mutex& mtx = conn.mutex();
		mtx.lock();
		ISMRMRD::Dataset dataset(filename.c_str(), groupname.c_str());
		mtx.unlock();
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
		for (i = images_.begin(); i != images_.end(); i++) {
			boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			ImageWrap& iw = *sptr_iw;
			int type = iw.type();
			void* ptr = iw.ptr_image();
			if (type == ISMRMRD::ISMRMRD_USHORT)
				writeImage(*(ISMRMRD::Image<unsigned short>*)ptr, dataset, mtx);
			else if (type == ISMRMRD::ISMRMRD_SHORT)
				writeImage(*(ISMRMRD::Image<short>*)ptr, dataset, mtx);
			else if (type == ISMRMRD::ISMRMRD_UINT)
				writeImage(*(ISMRMRD::Image<unsigned int>*)ptr, dataset, mtx);
			else if (type == ISMRMRD::ISMRMRD_INT)
				writeImage(*(ISMRMRD::Image<int>*)ptr, dataset, mtx);
			else if (type == ISMRMRD::ISMRMRD_FLOAT)
				writeImage(*(ISMRMRD::Image<float>*)ptr, dataset, mtx);
			else if (type == ISMRMRD::ISMRMRD_DOUBLE)
				writeImage(*(ISMRMRD::Image<double>*)ptr, dataset, mtx);
			else if (type == ISMRMRD::ISMRMRD_CXFLOAT)
				writeImage(*(ISMRMRD::Image< std::complex<float> >*)ptr, dataset, mtx);
			else if (type == ISMRMRD::ISMRMRD_CXDOUBLE)
				writeImage(*(ISMRMRD::Image< std::complex<double> >*)ptr, dataset, mtx);
		}
	}
	void getImageDimensions(int im_num, int* dim) {
		if (im_num < 0 || im_num >= images_.size())
			dim[0] = dim[1] = dim[2] = 0;
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
		int count = 0;
		for (i = images_.begin(); i != images_.end() && count < im_num; i++)
			count++;
		boost::shared_ptr<ImageWrap>& sptr_iw = *i;
		ImageWrap& iw = *sptr_iw;
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
	void getImageDataAsDoubleArray(int im_num, double* data) {
		if (im_num < 0 || im_num >= images_.size())
			return;
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
		int count = 0;
		for (i = images_.begin(); i != images_.end() && count < im_num; i++)
			count++;
		boost::shared_ptr<ImageWrap>& sptr_iw = *i;
		ImageWrap& iw = *sptr_iw;
		int type = iw.type();
		void* ptr = iw.ptr_image();
		if (type == ISMRMRD::ISMRMRD_USHORT)
			getImageData(*(ISMRMRD::Image<unsigned short>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_SHORT)
			getImageData(*(ISMRMRD::Image<short>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_UINT)
			getImageData(*(ISMRMRD::Image<unsigned int>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_INT)
			getImageData(*(ISMRMRD::Image<int>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_FLOAT)
			getImageData(*(ISMRMRD::Image<float>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_DOUBLE)
			getImageData(*(ISMRMRD::Image<double>*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_CXFLOAT)
			getImageData(*(ISMRMRD::Image< std::complex<float> >*)ptr, data);
		else if (type == ISMRMRD::ISMRMRD_CXDOUBLE)
			getImageData(*(ISMRMRD::Image< std::complex<double> >*)ptr, data);
	}

private:
	std::list<boost::shared_ptr<ImageWrap> > images_;
	template<typename T>
	void writeImage
		(ISMRMRD::Image<T>& im, ISMRMRD::Dataset& dataset, boost::mutex& mtx)
	{
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
	void getImageDim(ISMRMRD::Image<T>& im, int* dim) {
		dim[0] = im.getMatrixSizeX();
		dim[1] = im.getMatrixSizeY();
		dim[2] = im.getMatrixSizeZ();
	}
	template<typename T>
	void getImageData(ISMRMRD::Image<T>& im, double* data) {
		long long int n = im.getMatrixSizeX();
		n *= im.getMatrixSizeY();
		n *= im.getMatrixSizeZ();
		T* ptr = im.getDataPtr();
		for (long long int i = 0; i < n; i++)
			data[i] = std::fabs(ptr[i]);
	}
};


#endif
