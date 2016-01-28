#ifndef GADGETRON_EXTENSIONS
#define GADGETRON_EXTENSIONS

#include <cmath>

#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>

#include "gadgetron_client.h"
#include "gadget_lib.h"

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
	std::list<boost::shared_ptr<ImageWrap> >& operator()() {
		return images_;
	}
	const std::list<boost::shared_ptr<ImageWrap> >& operator()() const {
		return images_;
	}
	int size() const {
		return (int)images_.size();
	}
	void write(std::string filename, std::string groupname)
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
#ifdef MSVC
		std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#endif
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
#ifdef MSVC
		std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#endif
		int count = 0;
		for (i = images_.begin(); i != images_.end() && count < im_num; i++)
			count++;
		boost::shared_ptr<ImageWrap>& sptr_iw = *i;
		ImageWrap& iw = *sptr_iw;
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

#define ABS(X) ((X < 0) ? -X : X)

	template<typename T>
	void getImageData(ISMRMRD::Image<T>& im, double* data) {
		long long int n = im.getMatrixSizeX();
		n *= im.getMatrixSizeY();
		n *= im.getMatrixSizeZ();
		T* ptr = im.getDataPtr();
		for (long long int i = 0; i < n; i++)
			data[i] = ABS(ptr[i]);
	}
	template<typename T>
	void getImageUnsignedData(ISMRMRD::Image<T>& im, double* data) {
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

class GadgetHandle {
public:
	GadgetHandle(std::string id, boost::shared_ptr<aGadget> sptr_g) : 
		id_(id), sptr_g_(sptr_g) {}
	std::string id() const {
		return id_;
	}
	aGadget& gadget() {
		return *sptr_g_.get();
	}
	const aGadget& gadget() const {
		return *sptr_g_.get();
	}
private:
	std::string id_;
	boost::shared_ptr<aGadget> sptr_g_;
};

class GadgetChain {
public:
	void add_reader(std::string id, boost::shared_ptr<aGadget> sptr_g) {
			readers_.push_back(boost::shared_ptr<GadgetHandle>
				(new GadgetHandle(id, sptr_g)));
	}
	void add_writer(std::string id, boost::shared_ptr<aGadget> sptr_g) {
		writers_.push_back(boost::shared_ptr<GadgetHandle>
			(new GadgetHandle(id, sptr_g)));
	}
	void add_gadget(std::string id, boost::shared_ptr<aGadget> sptr_g) {
		gadgets_.push_back(boost::shared_ptr<GadgetHandle>
			(new GadgetHandle(id, sptr_g)));
	}
	std::string xml() const {
		std::string xml_script("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
		xml_script += "<gadgetronStreamConfiguration xsi:schemaLocation=";
		xml_script += "\"http://gadgetron.sf.net/gadgetron gadgetron.xsd\"\n";
        	xml_script += "xmlns=\"http://gadgetron.sf.net/gadgetron\"\n";
        	xml_script += "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n\n";

#ifdef MSVC
		std::list<boost::shared_ptr<GadgetHandle> >::const_iterator gh;
#else
		typename std::list<boost::shared_ptr<GadgetHandle> >::const_iterator gh;
#endif
		for (gh = readers_.begin(); gh != readers_.end(); gh++)
			xml_script += gh->get()->gadget().xml() + '\n';
		for (gh = writers_.begin(); gh != writers_.end(); gh++)
			xml_script += gh->get()->gadget().xml() + '\n';
		for (gh = gadgets_.begin(); gh != gadgets_.end(); gh++)
			xml_script += gh->get()->gadget().xml() + '\n';
		xml_script += "</gadgetronStreamConfiguration>\n";

		return xml_script;
	}
	//int nreaders() const {
	//	return readers_.size();
	//}
	//int nwriters() const {
	//	return writers_.size();
	//}
private:
	std::list<boost::shared_ptr<GadgetHandle> > readers_;
	std::list<boost::shared_ptr<GadgetHandle> > writers_;
	std::list<boost::shared_ptr<GadgetHandle> > gadgets_;
};

class MRIReconstruction : public GadgetChain {
public:
	MRIReconstruction() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdAcqMsgReader),
		writer_(new IsmrmrdImgMsgWriter),
		sptr_images_(new ImagesList) {
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		ImagesList& images = *sptr_images_;
		//ImagesList& images = *sptr_images_.get();
		con_().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(images())));
	}

	void process(ISMRMRD::Dataset& input) {

		std::string config = xml();
		//std::cout << config << std::endl;

		con_().connect(host_, port_);

		con_().send_gadgetron_configuration_script(config);

		input.readHeader(par_);
		con_().send_gadgetron_parameters(par_);

		boost::mutex& mtx = con_.mutex();
		uint32_t acquisitions = 0;
		{
			mtx.lock();
			acquisitions = input.getNumberOfAcquisitions();
			mtx.unlock();
		}

		//std::cout << acquisitions << " acquisitions" << std::endl;

		ISMRMRD::Acquisition acq_tmp;
		for (uint32_t i = 0; i < acquisitions; i++) {
			{
				boost::mutex::scoped_lock scoped_lock(mtx);
				input.readAcquisition(i, acq_tmp);
			}
			con_().send_ismrmrd_acquisition(acq_tmp);
		}

		con_().send_gadgetron_close();
		con_().wait();
	}

	boost::shared_ptr<ImagesList> get_output() {
		return sptr_images_;
	}

private:
	std::string host_;
	std::string port_;
	GTConnector con_;
	std::string par_;
	boost::shared_ptr<IsmrmrdAcqMsgReader> reader_;
	boost::shared_ptr<IsmrmrdImgMsgWriter> writer_;
	// for writing to a file
	//boost::shared_ptr<ISMRMRD::Dataset> sptr_images_ds_;
	boost::shared_ptr<ImagesList> sptr_images_;
};

class ImagesProcessor : public GadgetChain {
public:
	ImagesProcessor() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdImgMsgReader),
		writer_(new IsmrmrdImgMsgWriter),
		sptr_images_(new ImagesList) {
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		ImagesList& images = *sptr_images_;
		con_().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(images())));
	}

	void process(ImagesList& input) {
		std::string config = xml();
		//std::cout << config << std::endl;

		con_().connect(host_, port_);

		con_().send_gadgetron_configuration_script(config);

		std::list<boost::shared_ptr<ImageWrap> >& images = input();
#ifdef MSVC
		std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#else
		typename std::list<boost::shared_ptr<ImageWrap> >::iterator i;
#endif
		for (i = images.begin(); i != images.end(); i++) {
			boost::shared_ptr<ImageWrap>& sptr_iw = *i;
			ImageWrap& iw = *sptr_iw;
			con_().send_wrapped_image(iw);
		}

		con_().send_gadgetron_close();
		con_().wait();
	}

	boost::shared_ptr<ImagesList> get_output() {
		return sptr_images_;
	}

private:
	std::string host_;
	std::string port_;
	GTConnector con_;
	boost::shared_ptr<IsmrmrdImgMsgReader> reader_;
	boost::shared_ptr<IsmrmrdImgMsgWriter> writer_;
	boost::shared_ptr<ImagesList> sptr_images_;
};

#endif
