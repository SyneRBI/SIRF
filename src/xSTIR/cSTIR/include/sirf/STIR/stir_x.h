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

#ifndef EXTRA_STIR_TYPES
#define EXTRA_STIR_TYPES

#define WIN32_LEAN_AND_MEAN

/*!
\file
\ingroup STIR Extensions
\brief Specification file for extended STIR functionality classes.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#include <stdlib.h>

#include "sirf/STIR/stir_data_containers.h"

#define MIN_BIN_EFFICIENCY 1.0e-20f
//#define MIN_BIN_EFFICIENCY 1.0e-6f
#define DYNAMIC_CAST(T, X, Y) T& X = dynamic_cast<T&>(Y)

namespace sirf {

	/*!
\ingroup STIR Extensions
\brief Listmode-to-sinograms converter.

This class reads list mode data and produces corresponding *sinograms*,
i.e. histogrammed data in the format of PETAcquisitionData.

It has 2 main functions:
- process() can be used to read prompts and/or delayed coincidences to produce a single
PETAcquisitionData. 2 variables decide what done with 3 possible cases:
- `store_prompts`=`true`, `store_delayeds`=`false`: only prompts are stored
- `store_prompts`=`false`, `store_delayeds`=`true`: only delayeds are stored
- `store_prompts`=`true`, `store_delayeds`=`true`: prompts-delayeds are stored
Clearly, enabling the `store_delayeds` option only makes sense if the data was
acquired accordingly.
- estimate_randoms() can be used to get a relatively noiseless estimate of the
random coincidences.

Currently, the randoms are estimated from the delayed coincidences using the following
strategy:
1. singles (one per detector) are estimated using a Maximum Likelihood estimator
2. randoms-from-singles are computed per detector-pair via the usual product formula.
These are then added together for all detector pairs in a certain histogram-bin in the
data (accommodating for view mashing and axial compression).

The actual algorithm is described in

> D. Hogg, K. Thielemans, S. Mustafovic, and T. J. Spinks,
> "A study of bias for various iterative reconstruction methods in PET,"
> in 2002 IEEE Nuclear Science Symposium Conference Record, vol. 3. IEEE, Nov. 2002, pp. 1519-1523.
> [Online](http://dx.doi.org/10.1109/nssmic.2002.1239610).
*/

	class ListmodeToSinograms : public stir::LmToProjData {
	public:
		//! Constructor. 
		/*! Takes an optional text string argument with
			the name of a STIR parameter file defining the conversion options.
			If no argument is given, default settings apply except
			for the names of input raw data file, template file and
			output filename prefix, which must be set by the user by
			calling respective methods.

			By default, `store_prompts` is `true` and `store_delayeds` is `false`.
			*/
		//ListmodeToSinograms(const char* const par) : stir::LmToProjData(par) {}
		ListmodeToSinograms(const char* par) : stir::LmToProjData(par) {}
		ListmodeToSinograms() : stir::LmToProjData()
		{
			fan_size = -1;
			store_prompts = true;
			store_delayeds = false;
			delayed_increment = 0;
			num_iterations = 10;
			display_interval = 1;
			KL_interval = 1;
			save_interval = -1;
			//num_events_to_store = -1;
		}
		void set_input(std::string lm_file)
		{
			input_filename = lm_file;
		}
		//! Specifies the prefix for the output file(s), 
		/*! This will be appended by `_g1f1d0b0.hs`.
		*/
		void set_output(std::string proj_data_file)
		{
			output_filename_prefix = proj_data_file;
		}
		void set_template(std::string proj_data_file)
		{
			template_proj_data_name = proj_data_file;
		}
		void set_time_interval(double start, double stop)
		{
			std::pair<double, double> interval(start, stop);
			std::vector < std::pair<double, double> > intervals;
			intervals.push_back(interval);
			frame_defs = stir::TimeFrameDefinitions(intervals);
			do_time_frame = true;
		}
		int set_flag(const char* flag, bool value)
		{
			if (boost::iequals(flag, "store_prompts"))
				store_prompts = value;
			else if (boost::iequals(flag, "store_delayeds"))
				store_delayeds = value;
#if 0
			else if (boost::iequals(flag, "do_pre_normalisation"))
				do_pre_normalisation = value;
			else if (boost::iequals(flag, "do_time_frame"))
				do_time_frame = value;
#endif
			else if (boost::iequals(flag, "interactive"))
				interactive = value;
			else
				return -1;
			return 0;
		}
		bool get_store_prompts() const
		{
			return store_prompts;
		}
		bool get_store_delayeds() const
		{
			return store_delayeds;
		}
		bool set_up()
		{
			// always reset here, in case somebody set a new listmode or template file
			max_segment_num_to_process = -1;
			fan_size = -1;

			bool failed = post_processing();
			if (failed)
				return true;
			const int num_rings =
				lm_data_ptr->get_scanner_ptr()->get_num_rings();
			// code below is a copy of STIR for more generic cases
			if (max_segment_num_to_process == -1)
				max_segment_num_to_process = num_rings - 1;
			else
				max_segment_num_to_process =
				std::min(max_segment_num_to_process, num_rings - 1);
			const int max_fan_size =
				lm_data_ptr->get_scanner_ptr()->get_max_num_non_arccorrected_bins();
			if (fan_size == -1)
				fan_size = max_fan_size;
			else
				fan_size =
				std::min(fan_size, max_fan_size);
			half_fan_size = fan_size / 2;
			fan_size = 2 * half_fan_size + 1;

			return false;
		}
		stir::shared_ptr<PETAcquisitionData> get_output()
		{
			std::string filename = output_filename_prefix + "_f1g1d0b0.hs";
			return stir::shared_ptr<PETAcquisitionData>
				(new PETAcquisitionDataInFile(filename.c_str()));
		}

		int estimate_randoms()
		{
			compute_fan_sums_();
			int err = compute_singles_();
			if (err)
				return err;
			estimate_randoms_();
			return 0;
		}
		stir::shared_ptr<PETAcquisitionData> get_randoms_sptr()
		{
			return randoms_sptr;
		}
        /// Get the time at which the prompt rate exceeds a certain threshold.
        /// Returns -1 if not found.
        float get_time_at_which_prompt_rate_exceeds_threshold(const float threshold) const;

	protected:
		// variables for ML estimation of singles/randoms
		int fan_size;
		int half_fan_size;
		int max_ring_diff_for_fansums;
		int num_iterations;
		int display_interval;
		int KL_interval;
		int save_interval;
		stir::shared_ptr<std::vector<stir::Array<2, float> > > fan_sums_sptr;
		stir::shared_ptr<stir::DetectorEfficiencies> det_eff_sptr;
		stir::shared_ptr<PETAcquisitionData> randoms_sptr;
		void compute_fan_sums_(bool prompt_fansum = false);
		int compute_singles_();
		void estimate_randoms_();
		static unsigned long compute_num_bins_(const int num_rings,
			const int num_detectors_per_ring,
			const int max_ring_diff, const int half_fan_size);
	};

	/*!
	\ingroup STIR Extensions
	\brief Class for PET scanner detector efficiencies model.

	*/

	class PETAcquisitionSensitivityModel {
	public:
		PETAcquisitionSensitivityModel() {}
		// create from bin (detector pair) efficiencies sinograms
		PETAcquisitionSensitivityModel(PETAcquisitionData& ad);
		// create from ECAT8
		PETAcquisitionSensitivityModel(std::string filename);
		// chain two normalizations
		PETAcquisitionSensitivityModel
			(PETAcquisitionSensitivityModel& mod1, PETAcquisitionSensitivityModel& mod2)
		{
			norm_.reset(new stir::ChainedBinNormalisation(mod1.data(), mod2.data()));
		}

		stir::Succeeded set_up(const stir::shared_ptr<stir::ProjDataInfo>&);

		// multiply by bin efficiencies
		virtual void unnormalise(PETAcquisitionData& ad) const;
		// divide by bin efficiencies
		virtual void normalise(PETAcquisitionData& ad) const;
		// same as apply, but returns new data rather than changes old one
		stir::shared_ptr<PETAcquisitionData> forward(PETAcquisitionData& ad) const
		{
			stir::shared_ptr<PETAcquisitionData> sptr_ad = ad.new_acquisition_data();
			sptr_ad->fill(ad);
			this->unnormalise(*sptr_ad);
			return sptr_ad;
		}
		// same as undo, but returns new data rather than changes old one
		stir::shared_ptr<PETAcquisitionData> invert(PETAcquisitionData& ad) const
		{
			stir::shared_ptr<PETAcquisitionData> sptr_ad = ad.new_acquisition_data();
			sptr_ad->fill(ad);
			this->normalise(*sptr_ad);
			return sptr_ad;
		}

		stir::shared_ptr<stir::BinNormalisation> data()
		{
			return norm_;
			//return std::dynamic_pointer_cast<stir::BinNormalisation>(norm_);
		}

	protected:
		stir::shared_ptr<stir::BinNormalisation> norm_;
		//shared_ptr<stir::ChainedBinNormalisation> norm_;
	};

	
        /*!
	\ingroup STIR Extensions
	\brief A typedef to use SIRF terminology for DataProcessors

        \todo We should have a sirf::ImageDataProcessor which takes a sirf::ImageData, but that's too much work for now...
	*/
	typedef DataProcessor3DF ImageDataProcessor;

	/*!
	\ingroup STIR Extensions
	\brief Class for a PET acquisition model.

	PET acquisition model relates an image representation \e x to the
	acquisition data representation \e y as

	\f[ y = 1/n(G x + a) + b \f]

	where:
	<list>
	<item>
	\e G is the geometric (ray tracing) projector from the image voxels
	to the scanner's pairs of detectors (bins);
	</item>
	<item>
	\e a and \e b are otional additive and background terms representing
	the effects of noise and scattering; assumed to be 0 if not present;
	</item>
	<item>
	\e n is an optional bin normalization term representing the inverse of
	detector (bin) efficiencies; assumed to be 1 if not present.
	</item>
	</list>

	The computation of \e y for a given \e x by the above formula is
	referred to as forward projection, and the computation of

	\f[ z = G' m y \f]

	where \e G' is the transpose of \e G and \f$ m = 1/n \f$, is referred to as
	backward projection.

	There is a possibility to add an ImageDataProcessor to the acquisition model. Calling this
	\e P it extends the model to

	\f[ y = 1/n(G P x + a) + b \f]

	This can be used for instance to model resolution effects by first blurring the image.

	At present we use quick-fix implementation of forward projection for
	the computation of a subset of y. A more proper implementation will be done
	later via AcquisitionData subsets.
	*/

	class PETAcquisitionModel {
	public:
		void set_projectors(stir::shared_ptr<stir::ProjectorByBinPair> sptr_projectors)
		{
			sptr_projectors_ = sptr_projectors;
		}
		stir::shared_ptr<stir::ProjectorByBinPair> projectors_sptr()
		{
			return sptr_projectors_;
		}
		void set_additive_term(stir::shared_ptr<PETAcquisitionData> sptr)
		{
			sptr_add_ = sptr;
		}
		stir::shared_ptr<PETAcquisitionData> additive_term_sptr()
		{
			return sptr_add_;
		}
		void set_background_term(stir::shared_ptr<PETAcquisitionData> sptr)
		{
			sptr_background_ = sptr;
		}
		stir::shared_ptr<PETAcquisitionData> background_term_sptr()
		{
			return sptr_background_;
		}
		//void set_normalisation(shared_ptr<stir::BinNormalisation> sptr)
		//{
		//	sptr_normalisation_ = sptr;
		//}
		stir::shared_ptr<stir::BinNormalisation> normalisation_sptr()
		{
			if (sptr_asm_.get())
				return sptr_asm_->data();
			stir::shared_ptr<stir::BinNormalisation> sptr;
			return sptr;
			//return sptr_normalisation_;
		}
		//void set_bin_efficiency(shared_ptr<PETAcquisitionData> sptr_data);
		//void set_normalisation(shared_ptr<PETAcquisitionData> sptr_data)
		//{
		//	sptr_normalisation_.reset(new stir::BinNormalisationFromProjData(*sptr_data));
		//}
		void set_asm(stir::shared_ptr<PETAcquisitionSensitivityModel> sptr_asm)
		{
			//sptr_normalisation_ = sptr_asm->data();
			sptr_asm_ = sptr_asm;
		}

		//! sets data processor to use on the image before forward projection and after back projection
		/*! \warning This assumes that the data processor is its own adjoint.
		 */
		void set_image_data_processor(stir::shared_ptr<ImageDataProcessor> sptr_processor);
		void cancel_background_term()
		{
			sptr_background_.reset();
		}
		void cancel_additive_term()
		{
			sptr_add_.reset();
		}
		void cancel_normalisation()
		{
			sptr_asm_.reset();
			//sptr_normalisation_.reset();
		}

		virtual stir::Succeeded set_up(
			stir::shared_ptr<PETAcquisitionData> sptr_acq,
			stir::shared_ptr<STIRImageData> sptr_image);

		// computes and returns a subset of forward-projected data 
		stir::shared_ptr<PETAcquisitionData>
			forward(const STIRImageData& image,
			int subset_num = 0, int num_subsets = 1);
		// replaces a subset of acquisition data with forward-projected data
		void forward(PETAcquisitionData& acq_data, const STIRImageData& image,
			int subset_num, int num_subsets, bool zero = false);

		stir::shared_ptr<STIRImageData> backward(PETAcquisitionData& ad,
			int subset_num = 0, int num_subsets = 1);

	protected:
		stir::shared_ptr<stir::ProjectorByBinPair> sptr_projectors_;
		stir::shared_ptr<PETAcquisitionData> sptr_acq_template_;
		stir::shared_ptr<STIRImageData> sptr_image_template_;
		stir::shared_ptr<PETAcquisitionData> sptr_add_;
		stir::shared_ptr<PETAcquisitionData> sptr_background_;
		stir::shared_ptr<PETAcquisitionSensitivityModel> sptr_asm_;
		//shared_ptr<stir::BinNormalisation> sptr_normalisation_;
	};

	/*!
	\ingroup STIR Extensions
	\brief Ray tracing matrix implementation of the PET acquisition model.

	In this implementation \e x and \e y are essentially vectors and \e G
	a matrix. Each row of \e G corresponds to a line-of-response (LOR)
	between two detectors (there may be more than one line for each pair).
	The only non-zero elements of each row are those corresponding to
	voxels through which LOR passes, so the matrix is very sparse.
	Furthermore, owing to symmetries, many rows have the same values only
	in different order, and thus only one set of values needs to be computed
	and stored (see STIR documentation for details).
	*/

	class PETAcquisitionModelUsingMatrix : public PETAcquisitionModel {
	public:
		PETAcquisitionModelUsingMatrix()
		{
			this->sptr_projectors_.reset(new ProjectorPairUsingMatrix);
		}
		void set_matrix(stir::shared_ptr<stir::ProjMatrixByBin> sptr_matrix)
		{
			sptr_matrix_ = sptr_matrix;
			((ProjectorPairUsingMatrix*)this->sptr_projectors_.get())->
				set_proj_matrix_sptr(sptr_matrix);
		}
		stir::shared_ptr<stir::ProjMatrixByBin> matrix_sptr()
		{
			return ((ProjectorPairUsingMatrix*)this->sptr_projectors_.get())->
				get_proj_matrix_sptr();
		}
		virtual stir::Succeeded set_up(
			stir::shared_ptr<PETAcquisitionData> sptr_acq,
			stir::shared_ptr<STIRImageData> sptr_image)
		{
			if (!sptr_matrix_.get())
				return stir::Succeeded::no;
			return PETAcquisitionModel::set_up(sptr_acq, sptr_image);
		}

	private:
		stir::shared_ptr<stir::ProjMatrixByBin> sptr_matrix_;
	};

	typedef PETAcquisitionModel AcqMod3DF;
	typedef PETAcquisitionModelUsingMatrix AcqModUsingMatrix3DF;
	typedef stir::shared_ptr<AcqMod3DF> sptrAcqMod3DF;

#ifdef STIR_WITH_NIFTYPET_PROJECTOR
    /*!
    \ingroup STIR Extensions
    \brief NiftyPET implementation of the PET acquisition model.
    */

    class PETAcquisitionModelUsingNiftyPET : public PETAcquisitionModel {
    public:
        PETAcquisitionModelUsingNiftyPET()
        {
            _niftypet_projector_pair_sptr.reset(new ProjectorPairUsingNiftyPET);
            this->sptr_projectors_ = _niftypet_projector_pair_sptr;
        }
        void set_cuda_verbosity(const bool verbosity) const
        {
            _niftypet_projector_pair_sptr->set_verbosity(verbosity);
        }
    protected:
        stir::shared_ptr<ProjectorPairUsingNiftyPET> _niftypet_projector_pair_sptr;
    };
    typedef PETAcquisitionModelUsingNiftyPET AcqModUsingNiftyPET3DF;
#endif

	/*!
	\ingroup STIR Extensions
	\brief Attenuation model.

	*/

	class PETAttenuationModel : public PETAcquisitionSensitivityModel {
	public:
		PETAttenuationModel(STIRImageData& id, PETAcquisitionModel& am);
		// multiply by bin efficiencies
		virtual void unnormalise(PETAcquisitionData& ad) const;
		// divide by bin efficiencies
		virtual void normalise(PETAcquisitionData& ad) const;
	protected:
		stir::shared_ptr<stir::ForwardProjectorByBin> sptr_forw_projector_;
	};

	/*!
	\ingroup STIR Extensions
	\brief Accessor classes.

	Some methods of the STIR classes exposed to the user by SIRF are protected
	and hence cannot be called directly. The 'accessor' classes below bypass the
	protection by inheritance.
	*/

	class xSTIR_GeneralisedPrior3DF : public stir::GeneralisedPrior < Image3DF > {
	public:
		bool post_process() {
			return post_processing();
		}
	};

	class xSTIR_QuadraticPrior3DF : public stir::QuadraticPrior < float > {
	public:
		void only2D(int only) {
			only_2D = only != 0;
		}
	};

	class xSTIR_PLSPrior3DF : public stir::PLSPrior < float > {
	public:
		void only2D(int only) {
			only_2D = only != 0;
		}
	};

	class xSTIR_GeneralisedObjectiveFunction3DF :
		public stir::GeneralisedObjectiveFunction < Image3DF > {
	public:
		bool post_process() {
			return post_processing();
		}
	};

	//typedef xSTIR_GeneralisedObjectiveFunction3DF ObjectiveFunction3DF;

	class xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF :
		public stir::PoissonLogLikelihoodWithLinearModelForMeanAndProjData < Image3DF > {
	public:
		void set_input_file(const char* filename) {
			input_filename = filename;
		}
		void set_acquisition_data(stir::shared_ptr<PETAcquisitionData> sptr)
		{
			sptr_ad_ = sptr;
			set_proj_data_sptr(sptr->data());
		}
		void set_acquisition_model(stir::shared_ptr<AcqMod3DF> sptr)
		{
			sptr_am_ = sptr;
			AcqMod3DF& am = *sptr;
			set_projector_pair_sptr(am.projectors_sptr());
			if (am.additive_term_sptr().get())
				set_additive_proj_data_sptr(am.additive_term_sptr()->data());
			if (am.normalisation_sptr().get())
				set_normalisation_sptr(am.normalisation_sptr());
		}
		stir::shared_ptr<AcqMod3DF> acquisition_model_sptr()
		{
			return sptr_am_;
		}
	private:
		stir::shared_ptr<PETAcquisitionData> sptr_ad_;
		stir::shared_ptr<AcqMod3DF> sptr_am_;
	};

	typedef xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF
		PoissonLogLhLinModMeanProjData3DF;

	class xSTIR_IterativeReconstruction3DF :
		public stir::IterativeReconstruction < Image3DF > {
	public:
		bool post_process() {
			//std::cout << "in xSTIR_IterativeReconstruction3DF.post_process...\n";
			if (this->output_filename_prefix.length() < 1)
				this->set_output_filename_prefix("reconstructed_image");
			return post_processing();
		}
		stir::Succeeded setup(sptrImage3DF const& image) {
			//std::cout << "in xSTIR_IterativeReconstruction3DF.setup...\n";
			return set_up(image);
		}
		void update(Image3DF &image) {
			update_estimate(image);
			end_of_iteration_processing(image);
			subiteration_num++;
		}
		void update(STIRImageData& id)
		{
			update(id.data());
		}
		void update(stir::shared_ptr<STIRImageData> sptr_id)
		{
			update(*sptr_id);
		}
		int& subiteration() {
			return subiteration_num;
		}
		int subiteration() const {
			return subiteration_num;
		}
		void set_initial_estimate_file(const char* filename) {
			initial_data_filename = filename;
		}
	};

	class xSTIR_OSSPSReconstruction3DF : public stir::OSSPSReconstruction < Image3DF > {
	public:
		float& relaxation_parameter_value() {
			return relaxation_parameter;
		}
	};

	class xSTIR_FBP2DReconstruction : public stir::FBP2DReconstruction {
	public:
		xSTIR_FBP2DReconstruction()
		{
			_is_set_up = false;
		}
		void set_input(const PETAcquisitionData& acq)
		{
			set_input_data(acq.data());
		}
		void set_zoom(float v)
		{
                        set_zoom_xy(v);
		}
		void set_alpha_ramp(double alpha)
		{
			// does not work!
			//assert(alpha > 0 && alpha <= 1.0);
			if (!(alpha > 0 && alpha <= 1.0))
				throw LocalisedException
				("wrong ramp filter parameter alpha", __FILE__, __LINE__);
			alpha_ramp = alpha;
		}
		void set_frequency_cut_off(double fc)
		{
			if (!(fc > 0 && fc <= 0.5))
				throw LocalisedException
				("wrong frequency cut-off", __FILE__, __LINE__);
			fc_ramp = fc;
		}
		stir::Succeeded set_up(stir::shared_ptr<STIRImageData> sptr_id)
		{
			_sptr_image_data.reset(new STIRImageData(*sptr_id));
			_is_set_up = true;
			return stir::Succeeded::yes;
		}
		stir::Succeeded process()
		{
			if (!_is_set_up) {
				stir::shared_ptr<Image3DF> sptr_image(construct_target_image_ptr());
				_sptr_image_data.reset(new STIRImageData(sptr_image));
				return reconstruct(sptr_image);
			}
			else
				return reconstruct(_sptr_image_data->data_sptr());
		}
		stir::shared_ptr<STIRImageData> get_output()
		{
			return _sptr_image_data;
		}
	protected:
		bool _is_set_up;
		stir::shared_ptr<STIRImageData> _sptr_image_data;
	};

	class xSTIR_SeparableGaussianImageFilter : 
		public stir::SeparableGaussianImageFilter<float> {
	public:
		//stir::Succeeded set_up(const STIRImageData& id)
		//{
		//	return virtual_set_up(id.data());
		//}
		//void apply(STIRImageData& id)
		//{
		//	virtual_apply(id.data());
		//}
		void set_fwhms_xyz(char xyz, float f)
		{
			stir::BasicCoordinate<3, float> v = get_fwhms();
			switch (xyz) {
			case 'z':
				v[1] = f;
				break;
			case 'y':
				v[2] = f;
				break;
			case 'x':
				v[3] = f;
			}
			set_fwhms(v);
		}
		void set_max_kernel_sizes_xyz(char xyz, int s)
		{
			stir::BasicCoordinate<3, int> v = get_max_kernel_sizes();
			switch (xyz) {
			case 'z':
				v[1] = s;
				break;
			case 'y':
				v[2] = s;
				break;
			case 'x':
				v[3] = s;
			}
			set_max_kernel_sizes(v);
		}
	};
}

#endif
