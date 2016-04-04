#ifndef GADGETRON_EXTENSIONS
#define GADGETRON_EXTENSIONS

#include <cmath>

#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include "gadgetron_client.h"
#include "gadget_lib.h"
#include "ismrmrd_fftw.h"
#include "localised_exception.h"

#define N_TRIALS 5

class GTConnector {
public:
	GTConnector() {
		sptr_con_ = boost::shared_ptr<GadgetronClientConnector>
			(new GadgetronClientConnector);
	}
	GadgetronClientConnector& operator()() {
		return *sptr_con_.get();
	}
	boost::shared_ptr<GadgetronClientConnector> sptr() {
		return sptr_con_;
	}
private:
	boost::shared_ptr<GadgetronClientConnector> sptr_con_;
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
	virtual ~GadgetChain() 
	{
		//std::cout << "~GadgetChain called" << std::endl;
	}
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
	void set_endgadget(boost::shared_ptr<aGadget> sptr_g) {
		endgadget_ = sptr_g;
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
		xml_script += endgadget_->xml() + '\n';
		xml_script += "</gadgetronStreamConfiguration>\n";

		return xml_script;
	}
private:
	std::list<boost::shared_ptr<GadgetHandle> > readers_;
	std::list<boost::shared_ptr<GadgetHandle> > writers_;
	std::list<boost::shared_ptr<GadgetHandle> > gadgets_;
	boost::shared_ptr<aGadget> endgadget_;
};

class AcquisitionsProcessor : public GadgetChain {
public:
	AcquisitionsProcessor() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdAcqMsgReader),
		writer_(new IsmrmrdAcqMsgWriter)
	{
		sptr_acqs_.reset();
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		boost::shared_ptr<AcqFinishGadget> endgadget(new AcqFinishGadget);
		set_endgadget(endgadget);
	}
	virtual ~AcquisitionsProcessor() 
	{
		//std::cout << "~AcquisitionsProcessor called" << std::endl;
	}

	void process(AcquisitionsContainer& acquisitions) {

		std::string config = xml();
		//std::cout << config << std::endl;

		GTConnector conn;

		sptr_acqs_ = acquisitions.new_acquisitions_container();
		conn().register_reader(GADGET_MESSAGE_ISMRMRD_ACQUISITION,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientAcquisitionMessageCollector(sptr_acqs_)));

		for (int nt = 0; nt < N_TRIALS; nt++) {
			try {
				conn().connect(host_, port_);
				conn().send_gadgetron_configuration_script(config);

				conn().send_gadgetron_parameters(acquisitions.parameters());
				sptr_acqs_->copy_data(acquisitions);

				uint32_t nacq = 0;
				nacq = acquisitions.number();

				//std::cout << nacq << " acquisitions" << std::endl;

				ISMRMRD::Acquisition acq_tmp;
				for (uint32_t i = 0; i < nacq; i++) {
					acquisitions.get_acquisition(i, acq_tmp);
					conn().send_ismrmrd_acquisition(acq_tmp);
				}

				conn().send_gadgetron_close();
				conn().wait();

				break;
			}
			catch (...) {
				std::cout << "connection failed";
				if (nt < N_TRIALS - 1)
					std::cout << ", trying again...";
				std::cout << std::endl;
			}
		}
	}

	boost::shared_ptr<AcquisitionsContainer> get_output() {
		return sptr_acqs_;
	}

private:
	std::string host_;
	std::string port_;
	boost::shared_ptr<IsmrmrdAcqMsgReader> reader_;
	boost::shared_ptr<IsmrmrdAcqMsgWriter> writer_;
	boost::shared_ptr<AcquisitionsContainer> sptr_acqs_;
};

class ImagesReconstructor : public GadgetChain {
public:

	ImagesReconstructor() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdAcqMsgReader),
		writer_(new IsmrmrdImgMsgWriter)
	{
		sptr_images_.reset();
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		boost::shared_ptr<ImgFinishGadget> endgadget(new ImgFinishGadget);
		set_endgadget(endgadget);
	}

	void process(AcquisitionsContainer& acquisitions) {

		std::string config = xml();
		//std::cout << "config:\n" << config << std::endl;

		GTConnector conn;

		sptr_images_.reset(new ImagesList);
		conn().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(sptr_images_)));

		for (int nt = 0; nt < N_TRIALS; nt++) {
			try {
				conn().connect(host_, port_);
				conn().send_gadgetron_configuration_script(config);

				conn().send_gadgetron_parameters(acquisitions.parameters());

				uint32_t nacquisitions = 0;
				nacquisitions = acquisitions.number();

				//std::cout << nacquisitions << " acquisitions" << std::endl;

				ISMRMRD::Acquisition acq_tmp;
				for (uint32_t i = 0; i < nacquisitions; i++) {
					acquisitions.get_acquisition(i, acq_tmp);
					conn().send_ismrmrd_acquisition(acq_tmp);
				}

				conn().send_gadgetron_close();
				conn().wait();

				break;
			}
			catch (...) {
				std::cout << "connection failed";
				if (nt < N_TRIALS - 1)
					std::cout << ", trying again...";
				std::cout << std::endl;
			}
		}
	}

	boost::shared_ptr<ImagesContainer> get_output() {
		return sptr_images_;
	}

private:
	std::string host_;
	std::string port_;
	boost::shared_ptr<IsmrmrdAcqMsgReader> reader_;
	boost::shared_ptr<IsmrmrdImgMsgWriter> writer_;
	boost::shared_ptr<ImagesContainer> sptr_images_;
};

class ImagesProcessor : public GadgetChain {
public:
	ImagesProcessor() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdImgMsgReader),
		writer_(new IsmrmrdImgMsgWriter)
	{
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		boost::shared_ptr<ImgFinishGadget> endgadget(new ImgFinishGadget);
		set_endgadget(endgadget);
	}

	void process(ImagesContainer& images)
	{
		std::string config = xml();
		//std::cout << config << std::endl;

		GTConnector conn;

		sptr_images_ = images.new_images_container();
		conn().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(sptr_images_)));

		for (int nt = 0; nt < N_TRIALS; nt++) {
			try {
				conn().connect(host_, port_);
				conn().send_gadgetron_configuration_script(config);

				for (int i = 0; i < images.number(); i++) {
					ImageWrap& iw = images.image_wrap(i);
					conn().send_wrapped_image(iw);
				}

				conn().send_gadgetron_close();
				conn().wait();

				break;
			}
			catch (...) {
				std::cout << "connection failed";
				if (nt < N_TRIALS - 1)
					std::cout << ", trying again...";
				std::cout << std::endl;
			}
		}
	}

	boost::shared_ptr<ImagesContainer> get_output() {
		return sptr_images_;
	}

private:
	std::string host_;
	std::string port_;
	boost::shared_ptr<IsmrmrdImgMsgReader> reader_;
	boost::shared_ptr<IsmrmrdImgMsgWriter> writer_;
	boost::shared_ptr<ImagesContainer> sptr_images_;
};

class AcquisitionModel {
public:

	AcquisitionModel(
		boost::shared_ptr<AcquisitionsContainer> sptr_ac,
		boost::shared_ptr<ImagesContainer> sptr_ic
		) :
		sptr_acqs_(sptr_ac),
		sptr_imgs_(sptr_ic)
	{
		AcquisitionsContainer& ac = *sptr_ac;
		par_ = ac.parameters();
		ISMRMRD::deserialize(par_.c_str(), header_);
		ac.get_acquisition(0, acq_);
	}

	void setCSMs(boost::shared_ptr<CoilSensitivitiesContainer> sptr_csms)
	{
		sptr_csms_ = sptr_csms;
		//std::cout << sptr_csms_->items() << '\n';
	}

	void fwd(ImageWrap& iw, CoilSensitivityMap& csm, AcquisitionsContainer& ac)
	{
		int type = iw.type();
		void* ptr = iw.ptr_image();
		IMAGE_PROCESSING_SWITCH(type, fwd_, ptr, csm, ac);
	}

	void bwd(ImageWrap& iw, CoilSensitivityMap& csm, AcquisitionsContainer& ac, 
		int& off)
	{
		int type = iw.type();
		void* ptr = iw.ptr_image();
		IMAGE_PROCESSING_SWITCH(type, bwd_, ptr, csm, ac, off);
	}

	void fwd(ImagesContainer& ic, CoilSensitivitiesContainer& cc,
		AcquisitionsContainer& ac)
	{
		if (cc.items() < 1)
			throw LocalisedException
			("coil sensitivity maps not found", __FILE__, __LINE__);
		for (int i = 0; i < ic.number(); i++) {
			ImageWrap& iw = ic.image_wrap(i);
			CoilSensitivityMap& csm = cc(i%cc.items());
			fwd(iw, csm, ac);
		}
	}

	void bwd(ImagesContainer& ic, CoilSensitivitiesContainer& cc, 
		AcquisitionsContainer& ac)
	{
		if (cc.items() < 1)
			throw LocalisedException
			("coil sensitivity maps not found", __FILE__, __LINE__);
		ImageWrap iw(sptr_imgs_->image_wrap(0));
		//int dims[4];
		//iw.get_dim(dims);
		//for (int i = 0; i < ac.number() / dims[1]; i++) {
		for (int i = 0, a = 0; a < ac.number(); i++) {
			CoilSensitivityMap& csm = cc(i%cc.items());
			bwd(iw, csm, ac, a);
			ic.append(iw);
			//std::cout << i << ' ' << iw.norm() << std::endl;
		}
	}

	boost::shared_ptr<AcquisitionsContainer> fwd(ImagesContainer& ic)
	{
		if (!sptr_csms_.get() || sptr_csms_->items() < 1)
			throw LocalisedException
			("coil sensitivity maps not found", __FILE__, __LINE__);
		boost::shared_ptr<AcquisitionsContainer> sptr_acqs =
			sptr_acqs_->new_acquisitions_container();
		fwd(ic, *sptr_csms_, *sptr_acqs);
		return sptr_acqs;
	}

	boost::shared_ptr<ImagesContainer> bwd(AcquisitionsContainer& ac)
	{
		if (!sptr_csms_.get() || sptr_csms_->items() < 1)
			throw LocalisedException
			("coil sensitivity maps not found", __FILE__, __LINE__);
		boost::shared_ptr<ImagesContainer> sptr_imgs =
			sptr_imgs_->new_images_container();
		bwd(*sptr_imgs, *sptr_csms_, ac);
		return sptr_imgs;
	}

private:
	std::string par_;
	ISMRMRD::IsmrmrdHeader header_;
	ISMRMRD::Acquisition acq_;
	boost::shared_ptr<CoilSensitivitiesContainer> sptr_csms_;
	boost::shared_ptr<AcquisitionsContainer> sptr_acqs_;
	boost::shared_ptr<ImagesContainer> sptr_imgs_;

	float norm(ISMRMRD::NDArray<complex_float_t> arr)
	{
		float s = 0;
		complex_float_t* ia;
		for (ia = arr.begin(); ia != arr.end(); ia++)
			s += std::abs(std::conj(*ia) * (*ia));
		return sqrt(s);
	}

	template< typename T>
	void fwd_(ISMRMRD::Image<T>* ptr_img, CoilSensitivityMap& csm, 
		AcquisitionsContainer& ac)
	{
		ISMRMRD::Image<T>& img = *ptr_img;

		std::string par;
		ISMRMRD::IsmrmrdHeader header;
		par = ac.parameters();
		ISMRMRD::deserialize(par.c_str(), header);
		ISMRMRD::Encoding e = header.encoding[0];
		ISMRMRD::Acquisition acq(acq_);

		int readout = e.encodedSpace.matrixSize.x;
		unsigned int nx = e.reconSpace.matrixSize.x;
		unsigned int ny = e.reconSpace.matrixSize.y;
		unsigned int nc = acq.active_channels();

		//std::cout << nx << ' ' << ny << ' ' << nc << '\n';

		std::vector<size_t> dims;
		dims.push_back(readout); 
		dims.push_back(ny);
		dims.push_back(nc);

		ISMRMRD::NDArray<complex_float_t> cm(dims);
		memset(cm.getDataPtr(), 0, cm.getDataSize());

		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++) {
					uint16_t xout = x + (readout - nx) / 2;
					complex_float_t zi = (complex_float_t)img(x, y);
					complex_float_t zc = csm(x, y, 0, c);
					cm(xout, y, c) = zi * zc;
				}
			}
		}

		memset((void*)acq.getDataPtr(), 0, acq.getDataSize());

		fft2c(cm);

		for (size_t y = 0; y < ny; y++) {
			acq.clearAllFlags();
			if (y == 0)
				acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);
			if (y == ny - 1)
				acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
			acq.idx().kspace_encode_step_1 = y;
			acq.idx().repetition = 0;
			for (size_t c = 0; c < nc; c++) {
				for (size_t s = 0; s < readout; s++) {
					acq.data(s, c) = cm(s, y, c);
				}
			}
			ac.append_acquisition(acq);
		}
		ac.set_parameters(par);
		ac.write_data();

	}

	template< typename T>
	void bwd_(ISMRMRD::Image<T>* ptr_im, CoilSensitivityMap& csm, 
		AcquisitionsContainer& ac, int& off)
	{
		ISMRMRD::Image<T>& im = *ptr_im;

		std::string par;
		ISMRMRD::IsmrmrdHeader header;
		par = ac.parameters();
		ISMRMRD::deserialize(par.c_str(), header);
		ISMRMRD::Encoding e = header.encoding[0];
		ISMRMRD::Acquisition acq(acq_);

		int readout = e.encodedSpace.matrixSize.x;
		unsigned int nx = e.reconSpace.matrixSize.x;
		unsigned int ny = e.reconSpace.matrixSize.y;
		unsigned int nc = acq.active_channels();

		//std::cout << nx << ' ' << ny << ' ' << nc << '\n';

		std::vector<size_t> dims;
		dims.push_back(readout);
		dims.push_back(ny);
		dims.push_back(nc);

		ISMRMRD::NDArray<complex_float_t> cm(dims);
		for (size_t y = 0; y < ny; y++) {
			//ac.get_acquisition(y + ny*im_num, acq);
			ac.get_acquisition(y + off, acq);
			for (size_t c = 0; c < nc; c++) {
				for (size_t s = 0; s < readout; s++) {
					cm(s, y, c) = acq.data(s, c);
				}
			}
		}
		off += ny;
		ifft2c(cm);

		T* ptr = im.getDataPtr();
		T s;
		memset(ptr, 0, im.getDataSize());
		long long int i = 0;
		for (unsigned int c = 0; c < nc; c++) {
			i = 0;
			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++, i++) {
					uint16_t xout = x + (readout - nx) / 2;
					complex_float_t z = cm(xout, y, c);
					complex_float_t zc = csm(x, y, 0, c);
					xGadgetronUtilities::convert_complex(std::conj(zc) * z, s);
					ptr[i] += s;
				}
			}
		}

	}

};
#endif
