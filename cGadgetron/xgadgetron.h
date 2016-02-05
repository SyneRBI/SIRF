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
	boost::shared_ptr<aGadget> endgadget_;
};

class MRIReconstruction : public GadgetChain {
public:
	MRIReconstruction() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdAcqMsgReader),
		writer_(new IsmrmrdImgMsgWriter),
		sptr_images_(new ImagesList) 
	{
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		boost::shared_ptr<ImgFinishGadget> endgadget(new ImgFinishGadget);
		set_endgadget(endgadget);
		con_().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(sptr_images_)));
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

	void process(AcquisitionsContainer& acquisitions) {

		std::string config = xml();
		//std::cout << config << std::endl;

		con_().connect(host_, port_);

		con_().send_gadgetron_configuration_script(config);

		par_ = acquisitions.parameters();
		con_().send_gadgetron_parameters(par_);

		boost::mutex& mtx = con_.mutex();
		uint32_t nacquisitions = 0;
		{
			mtx.lock();
			nacquisitions = acquisitions.number();
			mtx.unlock();
		}

		//std::cout << nacquisitions << " acquisitions" << std::endl;

		ISMRMRD::Acquisition acq_tmp;
		for (uint32_t i = 0; i < nacquisitions; i++) {
			{
				boost::mutex::scoped_lock scoped_lock(mtx);
				acquisitions.getAcquisition(i, acq_tmp);
			}
			con_().send_ismrmrd_acquisition(acq_tmp);
		}

		con_().send_gadgetron_close();
		con_().wait();
	}

	boost::shared_ptr<ImagesContainer> get_output() {
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
		con_().register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(sptr_images_)));
	}

	void process(ImagesContainer& images)
	{
		std::string config = xml();
		//std::cout << config << std::endl;

		con_().connect(host_, port_);

		con_().send_gadgetron_configuration_script(config);

		for (int i = 0; i < images.number(); i++) {
			ImageWrap& iw = images.imageWrap(i);
			con_().send_wrapped_image(iw);
		}

		con_().send_gadgetron_close();
		con_().wait();
	}

	boost::shared_ptr<ImagesContainer> get_output() {
		return sptr_images_;
	}

private:
	std::string host_;
	std::string port_;
	GTConnector con_;
	boost::shared_ptr<IsmrmrdImgMsgReader> reader_;
	boost::shared_ptr<IsmrmrdImgMsgWriter> writer_;
	boost::shared_ptr<ImagesContainer> sptr_images_;
};

#endif
