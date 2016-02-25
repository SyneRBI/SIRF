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
	AcquisitionsProcessor(std::string filename) :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdAcqMsgReader),
		writer_(new IsmrmrdAcqMsgWriter),
		sptr_acqs_(new AcquisitionsFile(filename, true, true))
	{
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		boost::shared_ptr<AcqFinishGadget> endgadget(new AcqFinishGadget);
		set_endgadget(endgadget);
	}

	void process(AcquisitionsContainer& acquisitions) {

		std::string config = xml();
		//std::cout << config << std::endl;

		GTConnector conn;

		conn().register_reader(GADGET_MESSAGE_ISMRMRD_ACQUISITION,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientAcquisitionMessageCollector(sptr_acqs_)));

		conn().connect(host_, port_);
		conn().send_gadgetron_configuration_script(config);

		conn().send_gadgetron_parameters(acquisitions.parameters());
		sptr_acqs_->copyData(acquisitions);

		uint32_t nacq = 0;
		nacq = acquisitions.number();

		//std::cout << nacq << " acquisitions" << std::endl;

		ISMRMRD::Acquisition acq_tmp;
		for (uint32_t i = 0; i < nacq; i++) {
			acquisitions.getAcquisition(i, acq_tmp);
			conn().send_ismrmrd_acquisition(acq_tmp);
		}

		conn().send_gadgetron_close();
		conn().wait();
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

class ImageReconstructor : public GadgetChain {
public:

	ImageReconstructor() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdAcqMsgReader),
		writer_(new IsmrmrdImgMsgWriter),
		sptr_images_(new ImagesList) 
	{
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		boost::shared_ptr<ImgFinishGadget> endgadget(new ImgFinishGadget);
		set_endgadget(endgadget);
	}

	void process(AcquisitionsContainer& acquisitions) {

		std::string config = xml();
		//std::cout << "config:\n" << config << std::endl;

		GTConnector conn;

		conn().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(sptr_images_)));

		conn().connect(host_, port_);
		conn().send_gadgetron_configuration_script(config);

		conn().send_gadgetron_parameters(acquisitions.parameters());

		uint32_t nacquisitions = 0;
		nacquisitions = acquisitions.number();

		//std::cout << nacquisitions << " acquisitions" << std::endl;

		ISMRMRD::Acquisition acq_tmp;
		for (uint32_t i = 0; i < nacquisitions; i++) {
			acquisitions.getAcquisition(i, acq_tmp);
			conn().send_ismrmrd_acquisition(acq_tmp);
		}

		conn().send_gadgetron_close();
		conn().wait();
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
		writer_(new IsmrmrdImgMsgWriter),
		sptr_images_(new ImagesList) 
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

		conn().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(sptr_images_)));

		conn().connect(host_, port_);
		conn().send_gadgetron_configuration_script(config);

		for (int i = 0; i < images.number(); i++) {
			ImageWrap& iw = images.imageWrap(i);
			conn().send_wrapped_image(iw);
		}

		conn().send_gadgetron_close();
		conn().wait();
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
	AcquisitionModel(const AcquisitionsContainer& ac)
	{
		par_ = ac.parameters();
		coils_ = ac.coils();
		ISMRMRD::deserialize(par_.c_str(), header_);
	}

	void fwd(ImageWrap& iw, AcquisitionsContainer& ac)
	{
		int type = iw.type();
		void* ptr = iw.ptr_image();
		IMAGE_PROCESSING_SWITCH(type, fwd_, ptr, ac);
	}

private:
	std::string par_;
	ISMRMRD::IsmrmrdHeader header_;
	boost::shared_ptr<ISMRMRD::NDArray<complex_float_t> > coils_;

	template< typename T>
	void fwd_(ISMRMRD::Image<T>* ptr_im, AcquisitionsContainer& ac)
	{
		ISMRMRD::Image<T>& im = *ptr_im;
		ISMRMRD::Encoding e = header_.encoding[0];
		//ISMRMRD::AcquisitionSystemInformation sys = 
		//	header_.acquisitionSystemInformation;

		int readout = e.encodedSpace.matrixSize.x;
		unsigned int matrix_size = im.getMatrixSizeY();
		//unsigned int ncoils = sys.receiverChannels;
		const size_t *cdims = coils_->getDims();
		unsigned int ncoils = cdims[2];

		std::vector<size_t> dims;
		dims.push_back(readout); 
		dims.push_back(matrix_size);
		dims.push_back(ncoils);

		//boost::shared_ptr<NDArray<complex_float_t> > coils =
		//	generate_birdcage_sensititivies(matrix_size, ncoils, relative_radius_);
		ISMRMRD::NDArray<complex_float_t> coil_images(dims);
		memset(coil_images.getDataPtr(), 0, coil_images.getDataSize());

		T* ptr = im.getDataPtr();
		for (unsigned int c = 0; c < ncoils; c++) {
			long long int i = 0;
			for (unsigned int y = 0; y < matrix_size; y++) {
				for (unsigned int x = 0; x < matrix_size; x++, i++) {
					uint16_t xout = x + (readout - matrix_size) / 2;
					complex_float_t z = ptr[i];
					coil_images(xout, y, c) = z * (*coils_)(x, y, c);
				}
			}
		}

		ISMRMRD::Acquisition acq;
		acq.resize(readout, ncoils);
		memset((void*)acq.getDataPtr(), 0, acq.getDataSize());
		acq.available_channels() = ncoils;
		acq.center_sample() = (readout >> 1);

		ISMRMRD::NDArray<complex_float_t> cm = coil_images;
		fft2c(cm);
		for (size_t i = 0; i < matrix_size; i++) {
			acq.clearAllFlags();
			if (i == 0)
				acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);
			if (i == matrix_size - 1)
				acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
			acq.idx().kspace_encode_step_1 = i;
			acq.idx().repetition = 0;
			acq.sample_time_us() = 5.0;
			for (size_t c = 0; c < ncoils; c++) {
				for (size_t s = 0; s < readout; s++) {
					acq.data(s, c) = cm(s, i, c);
				}
			}
			ac.appendAcquisition(acq);
		}

	}
};
#endif
