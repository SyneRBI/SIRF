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

void
PETAcquisitionModel::set_bin_efficiency
(shared_ptr<PETAcquisitionData> sptr_data)
{
	////shared_ptr<ProjData> sptr(new ProjDataInMemory(*sptr_data));
	//shared_ptr<ProjData> sptr(new ProjDataInMemory(*sptr_data->data()));
	//inv_(sptr.get(), MIN_BIN_EFFICIENCY);
	//sptr_normalisation_.reset(new BinNormalisationFromProjData(sptr));

	shared_ptr<PETAcquisitionData>
		sptr_ad(sptr_data->new_acquisition_data());
	sptr_ad->inv(MIN_BIN_EFFICIENCY, *sptr_data);
	sptr_normalisation_.reset
		(new BinNormalisationFromProjData(sptr_ad->data()));
	sptr_normalisation_->set_up(sptr_ad->get_proj_data_info_sptr());

	//std::cout << sptr_normalisation_.get() << std::endl;
	//std::cout << sptr_normalisation_->is_trivial() << std::endl;
	//std::cout << normalisation_sptr()->is_trivial() << std::endl;
	//sptr_ad->clear_stream();
	sptr_norm_ = sptr_ad;
}

Succeeded 
PETAcquisitionModel::set_up(
	shared_ptr<PETAcquisitionData> sptr_acq,
	shared_ptr<Image3DF> sptr_image)
{
	Succeeded s = Succeeded::no;
	if (sptr_projectors_.get()) {
		s = sptr_projectors_->set_up
			(sptr_acq->get_proj_data_info_sptr(), sptr_image);
		sptr_acq_template_ = sptr_acq;
		sptr_image_template_ = sptr_image;
	}
	return s;
}

shared_ptr<PETAcquisitionData>
PETAcquisitionModel::forward(const Image3DF& image)
{
	shared_ptr<PETAcquisitionData> sptr_ad;
	sptr_ad = sptr_acq_template_->new_acquisition_data();

	shared_ptr<ProjData> sptr_fd = sptr_ad->data();

	sptr_projectors_->get_forward_projector_sptr()->forward_project
		(*sptr_fd, image);

	if (sptr_add_.get()) {
		std::cout << "additive term added...";
		sptr_ad->axpby(1.0, *sptr_ad, 1.0, *sptr_add_);
		//add_(sptr_fd, sptr_add_->data());
		std::cout << "ok\n";
	}
	else
		std::cout << "no additive term added\n";

	clear_stream();
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
		//add_(sptr_fd, sptr_background_->data());
		std::cout << "ok\n";
	}
	else
		std::cout << "no background term added\n";

	//sptr_ad->set_data(sptr_fd);
	//return sptr_fd;
	return sptr_ad;
}

shared_ptr<Image3DF> 
PETAcquisitionModel::backward(ProjData& ad)
{
	shared_ptr<Image3DF> sptr_im(sptr_image_template_->clone());
	sptr_im->fill(0.0);

	if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
		std::cout << "applying normalisation...";
		//ProjDataInMemory adc(ad);
		std::cout << "ok\n";
		sptr_normalisation_->undo(ad, 0, 1);
		//sptr_normalisation_->undo(adc, 0, 1);
		std::cout << "backprojecting...";
		sptr_projectors_->get_back_projector_sptr()->back_project(*sptr_im, ad);
		//(*sptr_im, adc);
		std::cout << "ok\n";
	}
	else
		sptr_projectors_->get_back_projector_sptr()->back_project
		(*sptr_im, ad);

	return sptr_im;
}

