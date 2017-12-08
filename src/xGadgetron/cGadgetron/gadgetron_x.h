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
\ingroup Gadgetron Extensions
\brief Specification file for extended Gadgetron functionality classes.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef GADGETRON_EXTENSIONS
#define GADGETRON_EXTENSIONS

#include <cmath>
#include <string>

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
#include "cgadgetron_shared_ptr.h"

#define N_TRIALS 5

/*!
\ingroup Gadgetron Extensions
\brief Shared pointer wrap-up for GadgetronClientConnector.

*/

using namespace gadgetron;

class GTConnector {
public:
	GTConnector() 
	{
		sptr_con_ = shared_ptr<GadgetronClientConnector>
			(new GadgetronClientConnector);
	}
	GadgetronClientConnector& operator()() 
	{
		return *sptr_con_.get();
	}
	shared_ptr<GadgetronClientConnector> sptr() 
	{
		return sptr_con_;
	}
private:
	shared_ptr<GadgetronClientConnector> sptr_con_;
};

/*!
\ingroup Gadgetron Extensions
\brief Shared pointer wrap-up for the abstract gadget class aGadget.

*/

class GadgetHandle {
public:
	GadgetHandle(std::string id, shared_ptr<aGadget> sptr_g) : 
		id_(id), sptr_g_(sptr_g) {}
	std::string id() const 
	{
		return id_;
	}
	aGadget& gadget() 
	{
		return *sptr_g_.get();
	}
	const aGadget& gadget() const 
	{
		return *sptr_g_.get();
	}
	shared_ptr<aGadget> gadget_sptr()
	{
		return sptr_g_;
	}
private:
	std::string id_;
	shared_ptr<aGadget> sptr_g_;
};

/*!
\ingroup Gadgetron Extensions
\brief Gadget chain class.

Gadgetron gadget chains have the following form:

-
reader gadget 
(receives data from the client and passes it down the chain)
-
first gadget 
(processes data received from the reader gadget and passes the result on)
-
second gadget 
(processes data received from the first gadget and passes the result on)
-
...
-
last gadget
(as above)
-
finishing gadget
(prepares the final result for sending to the client)
-
writer gadget
(sends the final result to the client)
*/

class GadgetChain { //: public anObject {
public:
	//GadgetChain()
	//{
	//	class_ = "GadgetChain";
	//}
	static const char* class_name()
	{
		return "GadgetChain";
	}
	// apparently caused crash in linux
	//virtual ~GadgetChain() {}
	// adds reader gadget
	void add_reader(std::string id, shared_ptr<aGadget> sptr_g) 
	{
			readers_.push_back(shared_ptr<GadgetHandle>
				(new GadgetHandle(id, sptr_g)));
	}
	// adds writer gadget
	void add_writer(std::string id, shared_ptr<aGadget> sptr_g) 
	{
		writers_.push_back(shared_ptr<GadgetHandle>
			(new GadgetHandle(id, sptr_g)));
	}
	// sdds finishig gadget
	void set_endgadget(shared_ptr<aGadget> sptr_g) 
	{
		endgadget_ = sptr_g;
	}
	// adds any other gadget
	void add_gadget(std::string id, shared_ptr<aGadget> sptr_g)
	{
		gadgets_.push_back(shared_ptr<GadgetHandle>
			(new GadgetHandle(id, sptr_g)));
	}
	shared_ptr<aGadget> gadget_sptr(std::string id);
	// returns string containing the definition of the chain in xml format
	std::string xml() const;
private:
	std::list<shared_ptr<GadgetHandle> > readers_;
	std::list<shared_ptr<GadgetHandle> > writers_;
	std::list<shared_ptr<GadgetHandle> > gadgets_;
	shared_ptr<aGadget> endgadget_;
};

/*!
\ingroup Gadgetron Extensions
\brief A particular type of Gadget chain that has AcquisitionData 
on input and output.

*/

class AcquisitionsProcessor : public GadgetChain {
public:
	AcquisitionsProcessor() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdAcqMsgReader),
		writer_(new IsmrmrdAcqMsgWriter)
	{
		//class_ = "AcquisitionsProcessor";
		sptr_acqs_.reset();
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		shared_ptr<AcquisitionFinishGadget> 
			endgadget(new AcquisitionFinishGadget);
		set_endgadget(endgadget);
	}
	// apparently caused crash in linux
	//virtual ~AcquisitionsProcessor() {}
	static const char* class_name()
	{
		return "AcquisitionsProcessor";
	}

	void process(MRAcquisitionData& acquisitions);
	shared_ptr<MRAcquisitionData> get_output() 
	{
		return sptr_acqs_;
	}

private:
	std::string host_;
	std::string port_;
	shared_ptr<IsmrmrdAcqMsgReader> reader_;
	shared_ptr<IsmrmrdAcqMsgWriter> writer_;
	shared_ptr<MRAcquisitionData> sptr_acqs_;
};

/*!
\ingroup Gadgetron Extensions
\brief A particular type of Gadget chain that has AcquisitionData on input 
and ImageData on output.

*/

class ImagesReconstructor : public GadgetChain {
public:

	ImagesReconstructor() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdAcqMsgReader),
		writer_(new IsmrmrdImgMsgWriter)
	{
		//class_ = "ImagesReconstructor";
		sptr_images_.reset();
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		shared_ptr<ImageFinishGadget> endgadget(new ImageFinishGadget);
		set_endgadget(endgadget);
	}
	static const char* class_name()
	{
		return "ImagesReconstructor";
	}

	void process(MRAcquisitionData& acquisitions);
	shared_ptr<MRImageData> get_output() 
	{
		return sptr_images_;
	}

private:
	std::string host_;
	std::string port_;
	shared_ptr<IsmrmrdAcqMsgReader> reader_;
	shared_ptr<IsmrmrdImgMsgWriter> writer_;
	shared_ptr<MRImageData> sptr_images_;
};

/*!
\ingroup Gadgetron Extensions
\brief A particular type of Gadget chain that has ImageData 
on input and output.

*/

class ImagesProcessor : public GadgetChain {
public:
	ImagesProcessor() :
		host_("localhost"), port_("9002"),
		reader_(new IsmrmrdImgMsgReader),
		writer_(new IsmrmrdImgMsgWriter)
	{
		//class_ = "ImagesProcessor";
		add_reader("reader", reader_);
		add_writer("writer", writer_);
		shared_ptr<ImageFinishGadget> endgadget(new ImageFinishGadget);
		set_endgadget(endgadget);
	}
	static const char* class_name()
	{
		return "ImagesProcessor";
	}

	void process(MRImageData& images);
	shared_ptr<MRImageData> get_output() 
	{
		return sptr_images_;
	}

private:
	std::string host_;
	std::string port_;
	shared_ptr<IsmrmrdImgMsgReader> reader_;
	shared_ptr<IsmrmrdImgMsgWriter> writer_;
	shared_ptr<MRImageData> sptr_images_;
};

/*!
\ingroup Gadgetron Extensions
\brief A class for MR acquisition modelling.

MR Acquisition model is a mathematical model that represents 
MR scanner by an operator (generally non-linear) that maps 
a mathematical representation of the scanned object \e x into 
a mathematical representation of the predicted acquisition 
data \e y to be produced by the scanner after scanning this
object:
\f[
  y = A(x).
\f]
The application of A is referred to as <em> (forward) projection </em>, 
and the application of the complex transpose of the Frechet 
derivative of A as \e backprojection.

In SIRF, \e x is represented by an \e ImageContainer object and 
\e y by an \e AcquisitionContainer object.
The application of A (projection) involves multiplication by
<em>coil sensitivity maps</em>, 2D Fourier transform applied
<em>xy</em>-slice-wise and gathering readouts based on the index
stored by property \e kspace_encode_step_1 of the method idx() of
the object of class ISMRMRD::Acquisition recorded by MRAcquisitionModel
constructor as a template.
*/

class MRAcquisitionModel {
public:

	/*
	The constructor records, by copying shared pointers, the two supplied
	arguments as templates, to be used for obtaining scanner and image
	discretisation data.
	*/
	MRAcquisitionModel(
		shared_ptr<MRAcquisitionData> sptr_ac,
		shared_ptr<MRImageData> sptr_ic
		) : sptr_acqs_(sptr_ac), sptr_imgs_(sptr_ic)
	{
	}

	// Records the coil sensitivities maps to be used. 
	void setCSMs(shared_ptr<CoilSensitivitiesContainer> sptr_csms)
	{
		sptr_csms_ = sptr_csms;
	}

	// Forward projects one image item (typically xy-slice) into
	// respective readouts, and appends them to the AcquisitionContainer
	// passed as the last argument.
	void fwd(ImageWrap& iw, CoilData& csm, MRAcquisitionData& ac, 
		unsigned int& off)
	{
		int type = iw.type();
		void* ptr = iw.ptr_image();
		IMAGE_PROCESSING_SWITCH(type, fwd_, ptr, csm, ac, off);
	}

	// Backprojects a set of readouts corresponding to one image item
	// (typically xy-slice).
	void bwd(ImageWrap& iw, CoilData& csm, MRAcquisitionData& ac, 
		unsigned int& off)
	{
		int type = iw.type();
		void* ptr = iw.ptr_image();
		IMAGE_PROCESSING_SWITCH(type, bwd_, ptr, csm, ac, off);
	}

	// Forward projects the whole ImageContainer using
	// coil sensitivity maps in the second argument.
	void fwd(MRImageData& ic, CoilSensitivitiesContainer& cc,
		MRAcquisitionData& ac);

	// Backprojects the whole AcquisitionContainer using
	// coil sensitivity maps in the second argument.
	void bwd(MRImageData& ic, CoilSensitivitiesContainer& cc,
		MRAcquisitionData& ac);

	// Forward projects the whole ImageContainer using
	// coil sensitivity maps referred to by sptr_csms_.
	shared_ptr<MRAcquisitionData> fwd(MRImageData& ic)
	{
		if (!sptr_csms_.get() || sptr_csms_->items() < 1)
			throw LocalisedException
			("coil sensitivity maps not found", __FILE__, __LINE__);
		shared_ptr<MRAcquisitionData> sptr_acqs = 
			sptr_acqs_->new_acquisitions_container();
		sptr_acqs->copy_acquisitions_info(*sptr_acqs_);
		fwd(ic, *sptr_csms_, *sptr_acqs);
		return sptr_acqs;
	}

	// Backprojects the whole AcquisitionContainer using
	// coil sensitivity maps referred to by sptr_csms_.
	shared_ptr<MRImageData> bwd(MRAcquisitionData& ac)
	{
		if (!sptr_csms_.get() || sptr_csms_->items() < 1)
			throw LocalisedException
			("coil sensitivity maps not found", __FILE__, __LINE__);
		shared_ptr<MRImageData> sptr_imgs =
			sptr_imgs_->new_images_container();
		bwd(*sptr_imgs, *sptr_csms_, ac);
		return sptr_imgs;
	}

private:
	std::string acqs_info_;
	shared_ptr<MRAcquisitionData> sptr_acqs_;
	shared_ptr<MRImageData> sptr_imgs_;
	shared_ptr<CoilSensitivitiesContainer> sptr_csms_;

	template< typename T>
	void fwd_(ISMRMRD::Image<T>* ptr_img, CoilData& csm,
		MRAcquisitionData& ac, unsigned int& off);
	template< typename T>
	void bwd_(ISMRMRD::Image<T>* ptr_im, CoilData& csm,
		MRAcquisitionData& ac, unsigned int& off);
};

#endif
