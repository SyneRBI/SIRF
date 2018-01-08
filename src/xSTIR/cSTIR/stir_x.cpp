/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2017 University College London.

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

#include "stir_types.h"
#include "stir_x.h"

using stir::shared_ptr;

PETAcquisitionSensitivityModel::
PETAcquisitionSensitivityModel(PETAcquisitionData& ad)
{
	shared_ptr<PETAcquisitionData>
		sptr_ad(ad.new_acquisition_data());
	sptr_ad->inv(MIN_BIN_EFFICIENCY, ad);
	shared_ptr<BinNormalisation> 
		sptr_n(new BinNormalisationFromProjData(sptr_ad->data()));
	shared_ptr<BinNormalisation> sptr_0;
	norm_.reset(new ChainedBinNormalisation(sptr_n, sptr_0));
}

PETAcquisitionSensitivityModel::
PETAcquisitionSensitivityModel(PETImageData& id)
{
	shared_ptr<BinNormalisationFromAttenuationImage>
		sptr_n(new BinNormalisationFromAttenuationImage(id.data_sptr()));
	shared_ptr<BinNormalisation> sptr_0;
	norm_.reset(new ChainedBinNormalisation(sptr_n, sptr_0));
}

PETAcquisitionSensitivityModel::
PETAcquisitionSensitivityModel(std::string filename)
{
	shared_ptr<BinNormalisationFromECAT8> 
		sptr_n(new BinNormalisationFromECAT8(filename));
	shared_ptr<BinNormalisation> sptr_0;
	norm_.reset(new ChainedBinNormalisation(sptr_n, sptr_0));
}

void
PETAcquisitionSensitivityModel::apply(PETAcquisitionData& ad)
{
	BinNormalisation* norm = norm_.get();
	norm->undo(*ad.data(), 0, 1);
}

void
PETAcquisitionSensitivityModel::undo(PETAcquisitionData& ad)
{
	BinNormalisation* norm = norm_.get();
	norm->apply(*ad.data(), 0, 1);
}

void
PETAcquisitionModel::set_bin_efficiency
(shared_ptr<PETAcquisitionData> sptr_data)
{
	shared_ptr<PETAcquisitionData>
		sptr_ad(sptr_data->new_acquisition_data());
	sptr_ad->inv(MIN_BIN_EFFICIENCY, *sptr_data);
	sptr_normalisation_.reset
		(new BinNormalisationFromProjData(sptr_ad->data()));
	sptr_normalisation_->set_up(sptr_ad->get_proj_data_info_sptr());

	//sptr_norm_ = sptr_ad;
}

Succeeded 
PETAcquisitionModel::set_up(
	shared_ptr<PETAcquisitionData> sptr_acq,
	shared_ptr<PETImageData> sptr_image)
{
	Succeeded s = Succeeded::no;
	if (sptr_projectors_.get()) {
		s = sptr_projectors_->set_up
			(sptr_acq->get_proj_data_info_sptr(), sptr_image->data_sptr());
		sptr_acq_template_ = sptr_acq;
		sptr_image_template_ = sptr_image;
	}
	return s;
}

shared_ptr<PETAcquisitionData>
PETAcquisitionModel::forward(const PETImageData& image)
{
	shared_ptr<PETAcquisitionData> sptr_ad;
	sptr_ad = sptr_acq_template_->new_acquisition_data();
	shared_ptr<ProjData> sptr_fd = sptr_ad->data();

	sptr_projectors_->get_forward_projector_sptr()->forward_project
		(*sptr_fd, image.data());

	if (sptr_add_.get()) {
		std::cout << "additive term added...";
		sptr_ad->axpby(1.0, *sptr_ad, 1.0, *sptr_add_);
		std::cout << "ok\n";
	}
	else
		std::cout << "no additive term added\n";

	//clear_stream();
	if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
		std::cout << "normalisation applied...";
		sptr_normalisation_->undo(*sptr_fd, 0, 1);
		//sptr_normalisation_->apply(*sptr_fd, 0, 1);
		std::cout << "ok\n";
	}
	else
		std::cout << "no normalisation applied\n";

	if (sptr_background_.get()) {
		std::cout << "background term added...";
		sptr_ad->axpby(1.0, *sptr_ad, 1.0, *sptr_background_);
		std::cout << "ok\n";
	}
	else
		std::cout << "no background term added\n";

	return sptr_ad;
}

shared_ptr<PETImageData> 
PETAcquisitionModel::backward(PETAcquisitionData& ad)
{
	shared_ptr<PETImageData> sptr_id;
	sptr_id = sptr_image_template_->new_image_data();
	shared_ptr<Image3DF> sptr_im = sptr_id->data_sptr();

	if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
		std::cout << "applying normalisation...";
		shared_ptr<PETAcquisitionData> sptr_ad(ad.new_acquisition_data());
		sptr_ad->fill(ad);
		sptr_normalisation_->undo(*sptr_ad->data(), 0, 1);
		std::cout << "ok\n";
		std::cout << "backprojecting...";
		sptr_projectors_->get_back_projector_sptr()->back_project(*sptr_im, *sptr_ad);
		std::cout << "ok\n";
	}
	else {
		std::cout << "backprojecting...";
		sptr_projectors_->get_back_projector_sptr()->back_project(*sptr_im, ad);
		std::cout << "ok\n";
	}

	return sptr_id;
}
