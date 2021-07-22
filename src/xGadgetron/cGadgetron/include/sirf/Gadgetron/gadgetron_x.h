/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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
\ingroup MR
\brief Specification file for extended Gadgetron functionality classes.

\author Evgueni Ovtchinnikov
\author SyneRBI
*/

#ifndef GADGETRON_EXTENSIONS
#define GADGETRON_EXTENSIONS

#define WIN32_LEAN_AND_MEAN

#include <cmath>
#include <list>
#include <string>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>

#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"

#include "sirf/common/JacobiCG.h"
#include "sirf/iUtilities/LocalisedException.h"
#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/gadget_lib.h"
#include "sirf/Gadgetron/ismrmrd_fftw.h"

#include "sirf/Gadgetron/FourierEncoding.h"

#if GADGETRON_TOOLBOXES_AVAILABLE
    #include "sirf/Gadgetron/NonCartesianEncoding.h"
#endif

#define N_TRIALS 5

namespace sirf {

  class GadgetronClientConnector;

	/*!
	\ingroup MR
	\brief Shared pointer wrap-up for GadgetronClientConnector.
	*/
	class GTConnector {
	public:
		GTConnector();
		GadgetronClientConnector& operator()();
		gadgetron::shared_ptr<GadgetronClientConnector> sptr();
	private:
		gadgetron::shared_ptr<GadgetronClientConnector> sptr_con_;
	};

	/*!
	\ingroup MR
	\brief Shared pointer wrap-up for the abstract gadget class aGadget.

	*/

	class GadgetHandle {
	public:
		GadgetHandle(std::string id, gadgetron::shared_ptr<aGadget> sptr_g) :
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
		gadgetron::shared_ptr<aGadget> gadget_sptr()
		{
			return sptr_g_;
		}
	private:
		std::string id_;
		gadgetron::shared_ptr<aGadget> sptr_g_;
	};

	/*!
	\ingroup MR
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
		GadgetChain() : host_("localhost"), port_("9002")
		{
			//class_ = "GadgetChain";
		}
		static const char* class_name()
		{
			return "GadgetChain";
		}
		void set_host(const std::string host)
		{
			host_ = host;
		}
		void set_port(const std::string port)
		{
			port_ = port;
		}
		// apparently caused crash in linux
		//virtual ~GadgetChain() {}
		// adds reader gadget
		void add_reader(std::string id, gadgetron::shared_ptr<aGadget> sptr_g)
		{
			readers_.push_back(gadgetron::shared_ptr<GadgetHandle>
				(new GadgetHandle(id, sptr_g)));
		}
		// adds writer gadget
		void add_writer(std::string id, gadgetron::shared_ptr<aGadget> sptr_g)
		{
			writers_.push_back(gadgetron::shared_ptr<GadgetHandle>
				(new GadgetHandle(id, sptr_g)));
		}
		// sdds finishig gadget
		void set_endgadget(gadgetron::shared_ptr<aGadget> sptr_g)
		{
			endgadget_ = sptr_g;
		}
		// adds any other gadget
		void add_gadget(std::string id, gadgetron::shared_ptr<aGadget> sptr_g)
		{
			gadgets_.push_back(gadgetron::shared_ptr<GadgetHandle>
				(new GadgetHandle(id, sptr_g)));
		}
		gadgetron::shared_ptr<aGadget> gadget_sptr(std::string id);
		// returns string containing the definition of the chain in xml format
		std::string xml() const;
	protected:
		std::string host_;
		std::string port_;
	private:
		std::list<gadgetron::shared_ptr<GadgetHandle> > readers_;
		std::list<gadgetron::shared_ptr<GadgetHandle> > writers_;
		std::list<gadgetron::shared_ptr<GadgetHandle> > gadgets_;
		gadgetron::shared_ptr<aGadget> endgadget_;
	};

	/*!
	\ingroup MR
	\brief A particular type of Gadget chain that has AcquisitionData
	on input and output.

	*/

	class AcquisitionsProcessor : public GadgetChain {
	public:
		AcquisitionsProcessor() :
			reader_(new IsmrmrdAcqMsgReader),
			writer_(new IsmrmrdAcqMsgWriter)
		{
			//class_ = "AcquisitionsProcessor";
			sptr_acqs_.reset();
			add_reader("reader", reader_);
			add_writer("writer", writer_);
		}
		// apparently caused crash in linux
		//virtual ~AcquisitionsProcessor() {}
		static const char* class_name()
		{
			return "AcquisitionsProcessor";
		}

		void process(MRAcquisitionData& acquisitions);
		gadgetron::shared_ptr<MRAcquisitionData> get_output()
		{
			if(!sptr_acqs_->sorted())
			 	sptr_acqs_->sort_by_time();

			return sptr_acqs_;
		}

	private:
		gadgetron::shared_ptr<IsmrmrdAcqMsgReader> reader_;
		gadgetron::shared_ptr<IsmrmrdAcqMsgWriter> writer_;
		gadgetron::shared_ptr<MRAcquisitionData> sptr_acqs_;
	};

	/*!
	\ingroup MR
	\brief A particular type of Gadget chain that has AcquisitionData on input
	and ImageData on output.

	*/

	class ImagesReconstructor : public GadgetChain {
	public:

		ImagesReconstructor() :
			reader_(new IsmrmrdAcqMsgReader),
			writer_(new IsmrmrdImgMsgWriter)
		{
			//class_ = "ImagesReconstructor";
			sptr_images_.reset();
			add_reader("reader", reader_);
			add_writer("writer", writer_);
			gadgetron::shared_ptr<ImageFinishGadget> endgadget(new ImageFinishGadget);
			set_endgadget(endgadget);
		}
		static const char* class_name()
		{
			return "ImagesReconstructor";
		}

		void process(MRAcquisitionData& acquisitions);
		gadgetron::shared_ptr<GadgetronImageData> get_output()
		{
			return sptr_images_;
		}

	private:
		gadgetron::shared_ptr<IsmrmrdAcqMsgReader> reader_;
		gadgetron::shared_ptr<IsmrmrdImgMsgWriter> writer_;
		gadgetron::shared_ptr<GadgetronImageData> sptr_images_;
	};

	/*!
	\ingroup MR
	\brief A particular type of Gadget chain that has ImageData
	on input and output.

	*/

	class ImagesProcessor : public GadgetChain {
	public:
		ImagesProcessor(bool dicom = false, std::string prefix = "image") :
			dicom_(dicom), prefix_(prefix),
			reader_(new IsmrmrdImgMsgReader)
		{
			//class_ = "ImagesProcessor";
			//gadgetron::shared_ptr<ImageFinishGadget> endgadget;
			gadgetron::shared_ptr<aGadget> endgadget;
			if (dicom) {
				writer_.reset(new DicomImageMessageWriter);
				endgadget.reset(new DicomFinishGadget);
			}
			else {
				writer_.reset(new IsmrmrdImgMsgWriter);
				endgadget.reset(new ImageFinishGadget);
			}
			add_reader("reader", reader_);
			add_writer("writer", writer_);
			set_endgadget(endgadget);
		}
		static const char* class_name()
		{
			return "ImagesProcessor";
		}

		void check_connection();
		void process(const GadgetronImageData& images);
		gadgetron::shared_ptr<GadgetronImageData> get_output()
		{
			return sptr_images_;
		}

	private:
		bool dicom_;
		std::string prefix_;
		gadgetron::shared_ptr<IsmrmrdImgMsgReader> reader_;
		gadgetron::shared_ptr<ImageMessageWriter> writer_;
//		gadgetron::shared_ptr<IsmrmrdImgMsgWriter> writer_;
		gadgetron::shared_ptr<GadgetronImageData> sptr_images_;
	};

	/*!
	\ingroup MR
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
		/*!
		\ingroup MR
		\brief Class for the product of backward and forward projectors of the MR acquisition model.

		For a given GadgetronImageData object x, computes B(F(x)), where F(x) is the forward projection of x,
		and B(y) is the backprojection of MRAcquisitionData object y.
		*/
		class BFOperator : public Operator<GadgetronImageData> {
		public:
			BFOperator(std::shared_ptr<MRAcquisitionModel> sptr_am) : sptr_am_(sptr_am) {}
			virtual std::shared_ptr<GadgetronImageData> 
				apply(GadgetronImageData& image_data)
			{
				std::shared_ptr<MRAcquisitionData> sptr_fwd =
					sptr_am_->fwd(image_data);
				std::shared_ptr<GadgetronImageData> sptr_bwd =
					sptr_am_->bwd(*sptr_fwd);
				return sptr_bwd;
			}
		private:
			std::shared_ptr<MRAcquisitionModel> sptr_am_;
		};

		MRAcquisitionModel() {}
		/*
		The constructor records, by copying shared pointers, the two supplied
		arguments as templates, to be used for obtaining scanner and image
		discretisation data.
		*/
		MRAcquisitionModel(
			gadgetron::shared_ptr<MRAcquisitionData> sptr_ac,
			gadgetron::shared_ptr<GadgetronImageData> sptr_ic
			)
		{
            this->set_up(sptr_ac, sptr_ic);
		}
		MRAcquisitionModel(
			gadgetron::shared_ptr<MRAcquisitionData> sptr_ac,
			gadgetron::shared_ptr<GadgetronImageData> sptr_ic,
			gadgetron::shared_ptr<CoilSensitivitiesVector> sptr_csms,
			std::string acqs_info
			) : sptr_csms_(sptr_csms), acqs_info_(acqs_info)
		{
			this->set_up(sptr_ac, sptr_ic);
		}
		
		float norm()
		{
			gadgetron::shared_ptr<MRAcquisitionModel> sptr_am
				(new MRAcquisitionModel(sptr_acqs_, sptr_imgs_, sptr_csms_, acqs_info_));

			BFOperator bf(sptr_am);
			JacobiCG<complex_float_t> jcg;
			jcg.set_num_iterations(2);
			gadgetron::unique_ptr<GadgetronImageData> sptr_id = sptr_imgs_->clone();
			GadgetronImageData& image_data = *sptr_id;
			image_data.fill(1.0);
			float lmd = jcg.largest(bf, image_data);
			return std::sqrt(lmd);
		}

		// make sure ic contains "true" images (and not e.g. G-factors)
		void check_data_role(const GadgetronImageData& ic);

		// Records the acquisition template to be used. 
		void set_acquisition_template
			(gadgetron::shared_ptr<MRAcquisitionData> sptr_ac)
		{
			sptr_acqs_ = sptr_ac;
		}
		// Records the image template to be used. 
		void set_image_template
			(gadgetron::shared_ptr<GadgetronImageData> sptr_ic)
		{
			check_data_role(*sptr_ic);
			sptr_imgs_ = sptr_ic;
		}
		// Returns shared pointer to the acquisitions template used. 
		gadgetron::shared_ptr<const MRAcquisitionData> acq_template_sptr() const
		{
			return sptr_acqs_;
		}
		// Returns shared pointer to the images template used. 
		gadgetron::shared_ptr<const GadgetronImageData> image_template_sptr() const
		{
			return sptr_imgs_;
		}
		// Records the coil sensitivities maps to be used. 
        void set_csm(gadgetron::shared_ptr<CoilSensitivitiesVector> sptr_csms)
		{
			sptr_csms_ = sptr_csms;
		}

        void set_encoder(gadgetron::shared_ptr<sirf::FourierEncoding> sptr_enc)
        {
            sptr_enc_ = sptr_enc;
        }

		// Records templates
		void set_up(gadgetron::shared_ptr<MRAcquisitionData> sptr_ac, 
			gadgetron::shared_ptr<GadgetronImageData> sptr_ic);
		
		// Forward projects the whole ImageContainer using
		// coil sensitivity maps in the second argument.
        void fwd(GadgetronImageData& ic, CoilSensitivitiesVector& cc,
			MRAcquisitionData& ac);

		// Backprojects the whole AcquisitionContainer using
		// coil sensitivity maps in the second argument.
        void bwd(GadgetronImageData& ic, const CoilSensitivitiesVector& cc,
            const MRAcquisitionData& ac);

		// Forward projects the whole ImageContainer using
		// coil sensitivity maps referred to by sptr_csms_.
		gadgetron::shared_ptr<MRAcquisitionData> fwd(GadgetronImageData& ic)
		{

            if (!sptr_acqs_.get())
				throw LocalisedException
				("Acquisition data template not set", __FILE__, __LINE__);
			if (!sptr_csms_.get() || sptr_csms_->items() < 1)
				throw LocalisedException
				("Coil sensitivity maps not found", __FILE__, __LINE__);
			check_data_role(ic);
            gadgetron::unique_ptr<MRAcquisitionData> uptr_acqs =
                sptr_acqs_->clone();

            fwd(ic, *sptr_csms_, *uptr_acqs);

            return std::shared_ptr<MRAcquisitionData>(std::move(uptr_acqs));
		}

		// Backprojects the whole AcquisitionContainer using
		// coil sensitivity maps referred to by sptr_csms_.
        gadgetron::shared_ptr<GadgetronImageData> bwd(const MRAcquisitionData& ac)
		{
			if (!sptr_imgs_.get())
				throw LocalisedException
				("image data template not set", __FILE__, __LINE__);
			if (!sptr_csms_.get() || sptr_csms_->items() < 1)
				throw LocalisedException
				("coil sensitivity maps not found", __FILE__, __LINE__);
			gadgetron::shared_ptr<GadgetronImageData> sptr_imgs =
				sptr_imgs_->new_images_container();
			bwd(*sptr_imgs, *sptr_csms_, ac);
			return sptr_imgs;
		}

	private:
		std::string acqs_info_;
		gadgetron::shared_ptr<MRAcquisitionData> sptr_acqs_;
		gadgetron::shared_ptr<GadgetronImageData> sptr_imgs_;
        gadgetron::shared_ptr<CoilSensitivitiesVector> sptr_csms_;
        gadgetron::shared_ptr<FourierEncoding> sptr_enc_;
	};

}

#endif
