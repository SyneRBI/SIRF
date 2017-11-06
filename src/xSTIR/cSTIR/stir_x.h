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

#include <stdlib.h>

#include "stir_data_containers.h"
#include "stir_types.h"

using stir::shared_ptr;

#define MIN_BIN_EFFICIENCY 1.0e-20f
//#define MIN_BIN_EFFICIENCY 1.0e-6f

class PETAcquisitionModel {
public:
	~PETAcquisitionModel()
	{
		sptr_projectors_.reset();
		sptr_acq_template_.reset();
		sptr_image_template_.reset();
		sptr_add_.reset();
		sptr_background_.reset();
		sptr_normalisation_.reset();
		sptr_norm_.reset();
	}
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
		shared_ptr<ProjData> sptr(new ProjDataInMemory(*sptr_data));
		sptr_normalisation_.reset(new BinNormalisationFromProjData(sptr));
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
	void clear_stream()
	{
		if (sptr_norm_.get())
			sptr_norm_->clear_stream();
	}
	void close_stream()
	{
		if (sptr_norm_.get())
			sptr_norm_->close_stream();
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
	shared_ptr<PETAcquisitionData> sptr_norm_;
};

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

class xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF :
	public PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> {
public:
	void set_input_file(const char* filename) {
		input_filename = filename;
	}
	void set_acquisition_data(shared_ptr<PETAcquisitionData> sptr_ad)
	{
		set_proj_data_sptr(sptr_ad->data());
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
	shared_ptr<AcqMod3DF> sptr_am_;
};

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

class xSTIR_OSSPSReconstruction3DF : public OSSPSReconstruction < Image3DF > {
public:
	float& relaxation_parameter_value() {
		return relaxation_parameter;
	}
};

#endif