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

/*!
\file
\ingroup STIR Extensions
\brief Specification file for extended STIR functionality classes.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#include <stdlib.h>

#include "stir_data_containers.h"
#include "stir_types.h"

using stir::shared_ptr;

#define MIN_BIN_EFFICIENCY 1.0e-20f
//#define MIN_BIN_EFFICIENCY 1.0e-6f

/*!
\ingroup STIR Extensions
\brief Listmode-to-sinograms converter.

*/

class ListmodeToSinograms : public LmToProjData {
public:
	//ListmodeToSinograms(const char* const par) : LmToProjData(par) {}
	ListmodeToSinograms(const char* par) : LmToProjData(par) {}
	ListmodeToSinograms() : LmToProjData()
	{
		fan_size = -1;
		delayed_increment = -1;
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
		frame_defs = TimeFrameDefinitions(intervals);
	}
	int set_flag(const char* flag, bool value)
	{
		if (boost::iequals(flag, "store_prompts"))
			store_prompts = value;
		else if (boost::iequals(flag, "store_delayeds"))
			store_delayeds = value;
		else if (boost::iequals(flag, "do_pre_normalisation"))
			do_pre_normalisation = value;
		else if (boost::iequals(flag, "do_time_frame"))
			do_time_frame = value;
		else if (boost::iequals(flag, "interactive"))
			interactive = value;
		else
			return -1;
		return 0;
	}
	bool set_up()
	{
		bool failed = post_processing();
		if (failed)
			return true;
		const int num_rings =
			lm_data_ptr->get_scanner_ptr()->get_num_rings();
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
		return false;
	}
	shared_ptr<PETAcquisitionData> get_output()
	{
		std::string filename = output_filename_prefix + "_f1g1d0b0.hs";
		return shared_ptr<PETAcquisitionData>
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
	shared_ptr<PETAcquisitionData> get_randoms_sptr()
	{
		return randoms_sptr;
	}

protected:
	int fan_size;
	int num_iterations;
	int display_interval;
	int KL_interval;
	int save_interval;
	shared_ptr<std::vector<Array<2, float> > > fan_sums_sptr;
	shared_ptr<DetectorEfficiencies> det_eff_sptr;
	shared_ptr<PETAcquisitionData> randoms_sptr;
	void compute_fan_sums_();
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
		norm_.reset(new ChainedBinNormalisation(mod1.data(), mod2.data()));
	}

	Succeeded set_up(const shared_ptr<ProjDataInfo>&);

	// multiply by bin efficiencies
	virtual void unnormalise(PETAcquisitionData& ad) const;
	// divide by bin efficiencies
	virtual void normalise(PETAcquisitionData& ad) const;
	// same as apply, but returns new data rather than changes old one
	shared_ptr<PETAcquisitionData> forward(PETAcquisitionData& ad) const
	{
		shared_ptr<PETAcquisitionData> sptr_ad = ad.new_acquisition_data();
		sptr_ad->fill(ad);
		this->unnormalise(*sptr_ad);
		return sptr_ad;
	}
	// same as undo, but returns new data rather than changes old one
	shared_ptr<PETAcquisitionData> invert(PETAcquisitionData& ad) const
	{
		shared_ptr<PETAcquisitionData> sptr_ad = ad.new_acquisition_data();
		sptr_ad->fill(ad);
		this->normalise(*sptr_ad);
		return sptr_ad;
	}

	shared_ptr<BinNormalisation> data()
	{
		return std::dynamic_pointer_cast<BinNormalisation>(norm_);
	}

protected:
	shared_ptr<ChainedBinNormalisation> norm_;
};

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
*/

class PETAcquisitionModel {
public:
	void set_projectors(shared_ptr<ProjectorByBinPair> sptr_projectors)
	{
		sptr_projectors_ = sptr_projectors;
	}
	shared_ptr<ProjectorByBinPair> projectors_sptr()
	{
		return sptr_projectors_;
	}
	void set_additive_term(shared_ptr<PETAcquisitionData> sptr)
	{
		sptr_add_ = sptr;
	}
	shared_ptr<PETAcquisitionData> additive_term_sptr()
	{
		return sptr_add_;
	}
	void set_background_term(shared_ptr<PETAcquisitionData> sptr)
	{
		sptr_background_ = sptr;
	}
	shared_ptr<PETAcquisitionData> background_term_sptr()
	{
		return sptr_background_;
	}
	void set_normalisation(shared_ptr<BinNormalisation> sptr)
	{
		sptr_normalisation_ = sptr;
	}
	shared_ptr<BinNormalisation> normalisation_sptr()
	{
		return sptr_normalisation_;
	}
	void set_bin_efficiency(shared_ptr<PETAcquisitionData> sptr_data);
	void set_normalisation(shared_ptr<PETAcquisitionData> sptr_data)
	{
		sptr_normalisation_.reset(new BinNormalisationFromProjData(*sptr_data));
	}
	void set_normalisation(shared_ptr<PETAcquisitionSensitivityModel> sptr_asm)
	{
		sptr_normalisation_ = sptr_asm->data();
	}

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
		sptr_normalisation_.reset();
	}

	virtual Succeeded set_up(
		shared_ptr<PETAcquisitionData> sptr_acq,
		shared_ptr<PETImageData> sptr_image);

	shared_ptr<PETAcquisitionData>
		forward(const PETImageData& image);

	shared_ptr<PETImageData> backward(PETAcquisitionData& ad);

protected:
	shared_ptr<ProjectorByBinPair> sptr_projectors_;
	shared_ptr<PETAcquisitionData> sptr_acq_template_;
	shared_ptr<PETImageData> sptr_image_template_;
	shared_ptr<PETAcquisitionData> sptr_add_;
	shared_ptr<PETAcquisitionData> sptr_background_;
	shared_ptr<BinNormalisation> sptr_normalisation_;
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
	void set_matrix(shared_ptr<ProjMatrixByBin> sptr_matrix)
	{
		sptr_matrix_ = sptr_matrix;
		((ProjectorPairUsingMatrix*)this->sptr_projectors_.get())->
			set_proj_matrix_sptr(sptr_matrix);
	}
	shared_ptr<ProjMatrixByBin> matrix_sptr()
	{
		return ((ProjectorPairUsingMatrix*)this->sptr_projectors_.get())->
			get_proj_matrix_sptr();
	}
	virtual Succeeded set_up(
		shared_ptr<PETAcquisitionData> sptr_acq,
		shared_ptr<PETImageData> sptr_image)
	{
		if (!sptr_matrix_.get())
			return Succeeded::no;
		return PETAcquisitionModel::set_up(sptr_acq, sptr_image);
	}

private:
	shared_ptr<ProjMatrixByBin> sptr_matrix_;
};

typedef PETAcquisitionModel AcqMod3DF;
typedef PETAcquisitionModelUsingMatrix AcqModUsingMatrix3DF;
typedef shared_ptr<AcqMod3DF> sptrAcqMod3DF;

/*!
\ingroup STIR Extensions
\brief Attenuation model.

*/

class PETAttenuationModel : public PETAcquisitionSensitivityModel {
public:
	PETAttenuationModel(PETImageData& id, PETAcquisitionModel& am);
	// multiply by bin efficiencies
	virtual void unnormalise(PETAcquisitionData& ad) const;
	// divide by bin efficiencies
	virtual void normalise(PETAcquisitionData& ad) const;
protected:
	shared_ptr<ProjectorByBinPair> sptr_projectors_;
};

/*!
\ingroup STIR Extensions
\brief Accessor classes.

Some methods of the STIR classes exposed to the user by SIRF are protected 
and hence cannot be called directly. The 'accessor' classes below bypass the 
protection by inheritance.
*/

class xSTIR_GeneralisedPrior3DF : public GeneralisedPrior < Image3DF > {
public:
	bool post_process() {
		return post_processing();
	}
};

class xSTIR_QuadraticPrior3DF : public QuadraticPrior < float > {
public:
	void only2D(int only) {
		only_2D = only != 0;
	}
};

class xSTIR_GeneralisedObjectiveFunction3DF :
	public GeneralisedObjectiveFunction<Image3DF> {
public:
	bool post_process() {
		return post_processing();
	}
};

//typedef xSTIR_GeneralisedObjectiveFunction3DF ObjectiveFunction3DF;

class xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF :
	public PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> {
public:
	void set_input_file(const char* filename) {
		input_filename = filename;
	}
	void set_acquisition_data(shared_ptr<PETAcquisitionData> sptr)
	{
		sptr_ad_ = sptr;
		set_proj_data_sptr(sptr->data());
	}
	void set_acquisition_model(shared_ptr<AcqMod3DF> sptr)
	{
		sptr_am_ = sptr;
		AcqMod3DF& am = *sptr;
		set_projector_pair_sptr(am.projectors_sptr());
		if (am.additive_term_sptr().get())
			set_additive_proj_data_sptr(am.additive_term_sptr()->data());
		if (am.normalisation_sptr().get())
			set_normalisation_sptr(am.normalisation_sptr());
	}
	shared_ptr<AcqMod3DF> acquisition_model_sptr()
	{
		return sptr_am_;
	}
private:
	shared_ptr<PETAcquisitionData> sptr_ad_;
	shared_ptr<AcqMod3DF> sptr_am_;
};

typedef xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF
PoissonLogLhLinModMeanProjData3DF;

class xSTIR_IterativeReconstruction3DF :
	public IterativeReconstruction<Image3DF> {
public:
	bool post_process() {
		if (this->output_filename_prefix.length() < 1)
			this->set_output_filename_prefix("reconstructed_image");
		return post_processing();
	}
	Succeeded setup(sptrImage3DF const& image) {
		return set_up(image);
	}
	void update(Image3DF &image) {
		update_estimate(image);
		end_of_iteration_processing(image);
		subiteration_num++;
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

class xSTIR_OSMAPOSLReconstruction3DF : 
	public OSMAPOSLReconstruction < Image3DF > {
public:
	Succeeded set_up(shared_ptr<PETImageData> sptr_id)
	{
		Succeeded s = Succeeded::no;
		xSTIR_IterativeReconstruction3DF* ptr_r =
			(xSTIR_IterativeReconstruction3DF*)this;
		if (!ptr_r->post_process()) {
			s = ptr_r->setup(sptr_id->data_sptr());
			ptr_r->subiteration() = ptr_r->get_start_subiteration_num();
		}
		return s;
	}
	void update(PETImageData& id)
	{
		((xSTIR_IterativeReconstruction3DF*)this)->update(id.data());
	}
	void update(shared_ptr<PETImageData> sptr_id)
	{
		update(*sptr_id);
	}
};

typedef xSTIR_OSMAPOSLReconstruction3DF OSMAPOSLReconstruction3DF;

class xSTIR_OSSPSReconstruction3DF : public OSSPSReconstruction < Image3DF > {
public:
	float& relaxation_parameter_value() {
		return relaxation_parameter;
	}
};

#endif