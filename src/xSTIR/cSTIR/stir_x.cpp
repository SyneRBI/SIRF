/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2020 University College London.

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

#include "stir/common.h"
#include "stir/config.h"
#include "stir/data/randoms_from_singles.h"
#include "stir/error.h"
#include "stir/IO/stir_ecat_common.h"
#include "stir/is_null_ptr.h"
#include "stir/multiply_crystal_factors.h"
#include "stir/Verbosity.h"

#include "sirf/STIR/stir_x.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

#if defined(HAVE_HDF5)
#include "stir/IO/GEHDF5Wrapper.h"
#include "stir/data/SinglesRatesFromGEHDF5.h"
#include "stir/recon_buildblock/BinNormalisationFromGEHDF5.h"
using namespace GE;
using namespace RDF_HDF5;
#endif

#ifdef STIR_USE_LISTMODEDATA
    typedef ListModeData LMD;
    typedef ListRecord LMR;
#else
    typedef CListModeData LMD;
    typedef CListRecord LMR;
#endif

float ListmodeToSinograms::get_time_at_which_num_prompts_exceeds_threshold(const unsigned long threshold) const
{
    if (input_filename.empty())
        throw std::runtime_error("ListmodeToSinograms::get_time_at_which_num_prompts_exceeds_threshold: Filename missing");

    shared_ptr<LMD> lm_data_ptr
      (read_from_file<LMD>(input_filename));

    shared_ptr <LMR> record_sptr = lm_data_ptr->get_empty_record_sptr();
    LMR& record = *record_sptr;

    double current_time = -1;
    unsigned long num_prompts = 0UL;

    /// Time resolution is 1s
    const double time_resolution = 1;

    while (true) {
        // no more events in file for some reason
        if (lm_data_ptr->get_next_record(record) == Succeeded::no)
            return -1.f;

        if (record.is_time()) {

            const double new_time = record.time().get_time_in_secs();
            // For the very first time
            if (current_time < 0) {
                current_time = new_time;
                num_prompts=0UL;
            }
            // Otherwise, increment the time
            else if (new_time >= current_time+time_resolution) {
                current_time += time_resolution;
                num_prompts=0UL;
            }
        }
        // If we found a prompt, increment!
        if (record.is_event() && record.event().is_prompt())
            ++num_prompts;
        // If the threshold is exceeded, return the time.
        if (num_prompts > threshold)
            return float(current_time);
    }
}

void
ListmodeToSinograms::compute_fan_sums_(bool prompt_fansum)
{
	//*********** get Scanner details
#if STIR_VERSION < 060000
        const auto& scanner = *lm_data_ptr->get_scanner_ptr();
#else
        const auto& scanner = lm_data_ptr->get_scanner();
#endif
	const auto num_rings = scanner.get_num_rings();
	const auto num_detectors_per_ring = scanner.get_num_detectors_per_ring();

	//*********** Finally, do the real work

	CPUTimer timer;
	timer.start();

	double time_of_last_stored_event = 0;
	long num_stored_events = 0;
	Array<2, float> data_fan_sums(IndexRange2D(num_rings, num_detectors_per_ring));
	fan_sums_sptr.reset(new std::vector<Array<2, float> >);

	// go to the beginning of the binary data
	lm_data_ptr->reset();

	max_ring_diff_for_fansums = 60;
	if (scanner != Scanner(Scanner::Siemens_mMR))
	{
		warning("This is not mMR data. Assuming all possible ring differences are in the listmode file");
		max_ring_diff_for_fansums = num_rings - 1;
	}
	unsigned int current_frame_num = 1;
	{
		// loop over all events in the listmode file
		shared_ptr<LMR> record_sptr =
			lm_data_ptr->get_empty_record_sptr();
		LMR& record = *record_sptr;

		bool first_event = true;

		double current_time = 0;
		while (true)
		{
			if (lm_data_ptr->get_next_record(record) == Succeeded::no)
			{
				// no more events in file for some reason
				std::cout << "processed frame " << current_frame_num << '\n';
				fan_sums_sptr->push_back(data_fan_sums);
				//write_fan_sums(data_fan_sums, current_frame_num);
				break; //get out of while loop
			}
			if (record.is_time())
			{
				const double new_time = record.time().get_time_in_secs();
				if (new_time >= frame_defs.get_end_time(current_frame_num) &&
					frame_defs.get_end_time(current_frame_num) > 
					frame_defs.get_start_time(current_frame_num))
				{
					while (current_frame_num <= frame_defs.get_num_frames() &&
						new_time >= frame_defs.get_end_time(current_frame_num))
					{
						//write_fan_sums(data_fan_sums, current_frame_num++);
						std::cout << "processed frame " << current_frame_num << '\n';
						fan_sums_sptr->push_back(data_fan_sums);
						current_frame_num++;
						data_fan_sums.fill(0);
					}
					if (current_frame_num > frame_defs.get_num_frames())
						break; // get out of while loop
				}
				current_time = new_time;
			}
			else if (record.is_event() &&
				frame_defs.get_start_time(current_frame_num) <= current_time)
			{
				// do a consistency check with dynamic_cast first
				if (first_event &&
					dynamic_cast<const CListEventCylindricalScannerWithDiscreteDetectors*>
					(&record.event()) == 0)
					error("Currently only works for scanners with discrete detectors.");
				first_event = false;

				// see if we increment or decrement the value in the sinogram
				const int event_increment =
                    (record.event().is_prompt() == prompt_fansum)
                    ? 1
                    : 0;

				if (event_increment == 0)
					continue;

				DetectionPositionPair<> det_pos;
				// because of above consistency check, we can use static_cast here 
				// (saving a bit of time)
				static_cast<const CListEventCylindricalScannerWithDiscreteDetectors&>
					(record.event()).get_detection_position(det_pos);
				const int ra = det_pos.pos1().axial_coord();
				const int rb = det_pos.pos2().axial_coord();
				const int a = det_pos.pos1().tangential_coord();
				const int b = det_pos.pos2().tangential_coord();
				if (abs(ra - rb) <= max_ring_diff_for_fansums)
				{
					const int det_num_diff =
						(a - b + 3 * num_detectors_per_ring / 2) % num_detectors_per_ring;
					if (det_num_diff <= fan_size / 2 ||
						det_num_diff >= num_detectors_per_ring - fan_size / 2)
					{
						data_fan_sums[ra][a] += event_increment;
						data_fan_sums[rb][b] += event_increment;
						num_stored_events += event_increment;
					}
					else
					{
					}
				}
				else
				{
				}

			} // end of spatial event processing
		} // end of while loop over all events

		time_of_last_stored_event =
			std::max(time_of_last_stored_event, current_time);
	}


	timer.stop();

	std::cerr << "Last stored event was recorded after time-tick at "
		<< time_of_last_stored_event << " secs\n";
	if (current_frame_num <= frame_defs.get_num_frames())
		std::cerr << "Early stop due to EOF. " << std::endl;
	std::cerr << "Total number of prompts/trues/delayed stored: "
		<< num_stored_events << std::endl;
	std::cerr << "\nThis took " << timer.value() << "s CPU time." << std::endl;

}

unsigned long
ListmodeToSinograms::compute_num_bins_(const int num_rings,
const int num_detectors_per_ring,
const int max_ring_diff, const int half_fan_size)
{
	unsigned long num = 0;
	for (int ra = 0; ra < num_rings; ++ra)
		for (int a = 0; a < num_detectors_per_ring; ++a)
		{
			for (int rb = std::max(ra - max_ring_diff, 0);
				rb <= std::min(ra + max_ring_diff, num_rings - 1); ++rb)
				for (int b = a + num_detectors_per_ring / 2 - half_fan_size;
					b <= a + num_detectors_per_ring / 2 + half_fan_size; ++b)
					++num;
		}
	return num;
}

int
ListmodeToSinograms::compute_singles_()
{
	const int do_display_interval = display_interval;
	const int do_KL_interval = KL_interval;
	const int do_save_interval =
		save_interval > 0 ? save_interval : num_iterations;

	int num_rings;
	int num_detectors_per_ring;
	int max_ring_diff = max_ring_diff_for_fansums;
	Array<2, float> data_fan_sums = (*fan_sums_sptr)[0];

	num_rings = data_fan_sums.get_length();
	ASSERT(num_rings > 0, "num_rings must be positive");
	ASSERT(data_fan_sums.get_min_index() == 0, "data_fan_sums.get_min_index() must be 0");
//	assert(num_rings > 0);
//	assert(data_fan_sums.get_min_index() == 0);
	num_detectors_per_ring = data_fan_sums[0].get_length();
	ASSERT(num_detectors_per_ring > 0, "num_detectors_per_ring must be positive");
//	assert(num_detectors_per_ring > 0);
	if (num_rings < max_ring_diff || num_detectors_per_ring < fan_size)
	{
		warning("fan sums matrix has sizes %dx%d, but this is "
			"too small compared to max_ring_diff (%d) and/or fan_size (%d)\n",
			num_rings, num_detectors_per_ring,
			max_ring_diff, fan_size);
		return EXIT_FAILURE;
	}

	CPUTimer timer;
	timer.start();

	det_eff_sptr.reset(new DetectorEfficiencies(
		IndexRange2D(num_rings, num_detectors_per_ring)));
	DetectorEfficiencies& efficiencies = *det_eff_sptr;
	//DetectorEfficiencies efficiencies(IndexRange2D(num_rings, num_detectors_per_ring));
	{
		float threshold_for_KL = data_fan_sums.find_max() / 100000.F;
		const int iter_num = 1;
		{
			if (iter_num == 1)
			{
				efficiencies.fill(sqrt(data_fan_sums.sum() /
					compute_num_bins_(num_rings, num_detectors_per_ring, max_ring_diff,
					half_fan_size)));
			}
			// efficiencies
			{
				for (int iter = 1; iter <= num_iterations; ++iter)
				{
					std::cout << "Starting iteration " << iter;
					iterate_efficiencies(efficiencies, data_fan_sums, max_ring_diff,
						half_fan_size);
					if (iter == num_iterations ||
						(do_KL_interval>0 && iter%do_KL_interval == 0))
					{
						Array<2, float> estimated_fan_sums(data_fan_sums.get_index_range());
						make_fan_sum_data(estimated_fan_sums, efficiencies, max_ring_diff,
							half_fan_size);
						std::cout << "\tKL " << KL(data_fan_sums, estimated_fan_sums,
							threshold_for_KL);
					}
					std::cout << std::endl;
				}
			} // end efficiencies

		}
	}
	timer.stop();
	std::cout << "CPU time " << timer.value() << " secs" << std::endl;
	return EXIT_SUCCESS;
}

int 
ListmodeToSinograms::estimate_randoms()
{
#if defined(HAVE_HDF5)
	std::cout << "estimate_randoms: trying GEHDF5...\n";
	try {
		if (GEHDF5Wrapper::check_GE_signature(input_filename)) {
			SinglesRatesFromGEHDF5  singles;
			singles.read_from_file(input_filename);
			GEHDF5Wrapper input_file(input_filename);
#if STIR_VERSION < 060000
			float coincidence_time_window = input_file.get_coincidence_time_window();
			ProjData& proj_data = *randoms_sptr->data();
			randoms_from_singles(proj_data, singles, coincidence_time_window);
#else
                        randoms_from_singles(*randoms_sptr->data(), singles);
#endif
                        return 0;
		}
	}
	catch (...) {
	}
	std::cout << "not a GE HDF5 file. Using ML estimate from delayeds\n";
#endif
	compute_fan_sums_();
	int err = compute_singles_();
	if (err)
		return err;
	ProjData& proj_data = *randoms_sptr->data();
	DetectorEfficiencies& efficiencies = *det_eff_sptr;
	multiply_crystal_factors(proj_data, efficiencies, 1.0f);
	return 0;
}

PETAcquisitionSensitivityModel::
PETAcquisitionSensitivityModel(STIRAcquisitionData& ad)
{
        std::shared_ptr<STIRAcquisitionData>
		sptr_ad(ad.new_acquisition_data());
	sptr_ad->inv(MIN_BIN_EFFICIENCY, ad);
        stir::shared_ptr<BinNormalisation> 
		sptr_n(new BinNormalisationFromProjData(sptr_ad->data()));
	//shared_ptr<BinNormalisation> sptr_0;
	//norm_.reset(new ChainedBinNormalisation(sptr_n, sptr_0));
	norm_ = sptr_n;
	//norm_ = shared_ptr<BinNormalisation>
	//	(new BinNormalisationFromProjData(sptr_ad->data()));
}

PETAcquisitionSensitivityModel::
PETAcquisitionSensitivityModel(std::string filename)
{
#if defined(HAVE_HDF5)
	std::cout << "trying GEHDF5...\n";
	try {
		if (GEHDF5Wrapper::check_GE_signature(filename)) {
			shared_ptr<BinNormalisation>
				sptr_n(new BinNormalisationFromGEHDF5(filename));
			norm_ = sptr_n;
			std::cout << "created bin normalisation from GE HDF5 file\n";
			return;
		}
	}
	catch (...) {
	}
	std::cout << "not a GE HDF5 file\n";
#endif
	shared_ptr<BinNormalisation>
		sptr_n(new BinNormalisationFromECAT8(filename));
	norm_ = sptr_n;
}

void
PETAcquisitionSensitivityModel::set_up(const shared_ptr<const ExamInfo>& sptr_ei,
	const shared_ptr<ProjDataInfo>& sptr_pdi)
{
#if STIR_VERSION < 050000
	stir::Succeeded s = norm_->set_up(sptr_pdi);
#else
	stir::Succeeded s = norm_->set_up(sptr_ei, sptr_pdi);
#endif
	if (s != stir::Succeeded::yes)
		THROW("stir::BinNormalisation setup failed");
}

void
PETAcquisitionSensitivityModel::unnormalise(STIRAcquisitionData& ad) const
{
	BinNormalisation* norm = norm_.get();
#if STIR_VERSION < 050000
	norm->undo(*ad.data(), 0, 1);
#else
	norm->undo(*ad.data());
#endif
}

void
PETAcquisitionSensitivityModel::normalise(STIRAcquisitionData& ad) const
{
	BinNormalisation* norm = norm_.get();
#if STIR_VERSION < 050000
	norm->apply(*ad.data(), 0, 1);
#else
	norm->apply(*ad.data());
#endif
}

PETAttenuationModel::PETAttenuationModel
(STIRImageData& id, PETAcquisitionModel& am)
{
	sptr_forw_projector_ = am.projectors_sptr()->get_forward_projector_sptr();
        if (is_null_ptr(sptr_forw_projector_))
          error("PETAttenuationModel: Forward projector not set correctly. Something wrong.");
	shared_ptr<BinNormalisation>
		sptr_n(new BinNormalisationFromAttenuationImage
		(id.data_sptr(), sptr_forw_projector_));
	norm_ = sptr_n;
}

void
PETAttenuationModel::unnormalise(STIRAcquisitionData& ad) const
{
	//std::cout << "in PETAttenuationModel::unnormalise\n";
	BinNormalisation* norm = norm_.get();
        stir::shared_ptr<DataSymmetriesForViewSegmentNumbers>
		symmetries_sptr(sptr_forw_projector_->get_symmetries_used()->clone());
#if STIR_VERSION < 050000
	norm->undo(*ad.data(), 0, 1, symmetries_sptr);
#else
	norm->undo(*ad.data(), symmetries_sptr);
#endif
}

void
PETAttenuationModel::normalise(STIRAcquisitionData& ad) const
{
	BinNormalisation* norm = norm_.get();
        stir::shared_ptr<DataSymmetriesForViewSegmentNumbers>
		symmetries_sptr(sptr_forw_projector_->get_symmetries_used()->clone());
#if STIR_VERSION < 050000
	norm->apply(*ad.data(), 0, 1, symmetries_sptr);
#else
	norm->apply(*ad.data(), symmetries_sptr);
#endif
}

//void
//PETAcquisitionModel::set_bin_efficiency
//(shared_ptr<STIRAcquisitionData> sptr_data)
//{
//	std::shared_ptr<STIRAcquisitionData>
//		sptr_ad(sptr_data->new_acquisition_data());
//	sptr_ad->inv(MIN_BIN_EFFICIENCY, *sptr_data);
//	sptr_normalisation_.reset
//		(new BinNormalisationFromProjData(sptr_ad->data()));
//	sptr_normalisation_->set_up(sptr_ad->get_proj_data_info_sptr());
//}

void
PETAcquisitionModel::set_up(
                            std::shared_ptr<STIRAcquisitionData> sptr_acq,
                            std::shared_ptr<STIRImageData> sptr_image)
{
	Succeeded s = Succeeded::no;
	if (sptr_projectors_.get()) {
		s = sptr_projectors_->set_up
			(sptr_acq->get_proj_data_info_sptr()->create_shared_clone(),
				sptr_image->data_sptr());
		sptr_acq_template_ = sptr_acq;
		sptr_image_template_ = sptr_image;
	}
	if (s == Succeeded(Succeeded::yes)) {
		if (sptr_asm_ && sptr_asm_->data())
			sptr_asm_->set_up(sptr_acq->get_exam_info_sptr(),
				sptr_acq->get_proj_data_info_sptr()->create_shared_clone());
	}
	else
		THROW("stir::ProjectorByBinPair setup failed");
}

void 
PETAcquisitionModel::set_image_data_processor(stir::shared_ptr<ImageDataProcessor> sptr_processor)
{
	if (!sptr_projectors_)
		throw std::runtime_error("projectors need to be set before calling set_image_data_processor");

	sptr_projectors_->get_forward_projector_sptr()->set_pre_data_processor(sptr_processor);
	sptr_projectors_->get_back_projector_sptr()->set_post_data_processor(sptr_processor);
}

void 
PETAcquisitionModel::forward(STIRAcquisitionData& ad, const STIRImageData& image,
	int subset_num, int num_subsets, bool zero, bool do_linear_only) const
{
        stir::shared_ptr<ProjData> sptr_fd = ad.data();
	sptr_projectors_->get_forward_projector_sptr()->forward_project
		(*sptr_fd, image.data(), subset_num, num_subsets, zero);

	float one = 1.0;

	if (sptr_add_.get() && !do_linear_only) {
		if (stir::Verbosity::get() > 1) std::cout << "additive term added...";
		ad.axpby(&one, ad, &one, *sptr_add_);
		//ad.axpby(1.0, ad, 1.0, *sptr_add_);
		if (stir::Verbosity::get() > 1) std::cout << "ok\n";
	}
	else
		if (stir::Verbosity::get() > 1) std::cout << "no additive term added\n";

	PETAcquisitionSensitivityModel* sm = sptr_asm_.get();
	if (sm && sm->data() && !sm->data()->is_trivial()) {
		if (stir::Verbosity::get() > 1) std::cout << "applying unnormalisation...";
		sptr_asm_->unnormalise(ad);
		if (stir::Verbosity::get() > 1) std::cout << "ok\n";
	}
	else
		if (stir::Verbosity::get() > 1) std::cout << "no unnormalisation applied\n";

	if (sptr_background_.get() && !do_linear_only) {
		if (stir::Verbosity::get() > 1) std::cout << "background term added...";
		ad.axpby(&one, ad, &one, *sptr_background_);
		//ad.axpby(1.0, ad, 1.0, *sptr_background_);
		if (stir::Verbosity::get() > 1) std::cout << "ok\n";
	}
	else
		if (stir::Verbosity::get() > 1) std::cout << "no background term added\n";
}

std::shared_ptr<STIRAcquisitionData>
PETAcquisitionModel::forward(const STIRImageData& image, 
	int subset_num, int num_subsets, bool do_linear_only) const
{
	if (!sptr_acq_template_.get())
		THROW("Fatal error in PETAcquisitionModel::forward: acquisition template not set");
        std::shared_ptr<STIRAcquisitionData> sptr_ad =
          sptr_acq_template_->new_acquisition_data();
        stir::shared_ptr<ProjData> sptr_fd = sptr_ad->data();
	forward(*sptr_ad, image, subset_num, num_subsets, num_subsets > 1, do_linear_only);
	return sptr_ad;
}

std::shared_ptr<STIRImageData> 
PETAcquisitionModel::backward(const STIRAcquisitionData& ad,
	int subset_num, int num_subsets) const
{
	if (!sptr_image_template_.get())
		THROW("Fatal error in PETAcquisitionModel::backward: image template not set");
        std::shared_ptr<STIRImageData> sptr_id;
	sptr_id = sptr_image_template_->new_image_data();
	backward(*sptr_id, ad, subset_num, num_subsets);
	return sptr_id;
}

void
PETAcquisitionModel::backward(STIRImageData& id, const STIRAcquisitionData& ad,
	int subset_num, int num_subsets) const
{
        stir::shared_ptr<Image3DF> sptr_im = id.data_sptr();

	PETAcquisitionSensitivityModel* sm = sptr_asm_.get();
	if (sm && sm->data() && !sm->data()->is_trivial()) {
		if (stir::Verbosity::get() > 1) std::cout << "applying unnormalisation...";
                std::shared_ptr<STIRAcquisitionData> sptr_ad(ad.new_acquisition_data());
		sptr_ad->fill(ad);
		sptr_asm_->unnormalise(*sptr_ad);
		//sptr_normalisation_->undo(*sptr_ad->data(), 0, 1);
		if (stir::Verbosity::get() > 1) std::cout << "ok\n";
		if (stir::Verbosity::get() > 1) std::cout << "backprojecting...";
		sptr_projectors_->get_back_projector_sptr()->back_project
			(*sptr_im, *sptr_ad, subset_num, num_subsets);
		if (stir::Verbosity::get() > 1) std::cout << "ok\n";
	}
	else {
		if (stir::Verbosity::get() > 1) std::cout << "backprojecting...";
		sptr_projectors_->get_back_projector_sptr()->back_project
			(*sptr_im, ad, subset_num, num_subsets);
		if (stir::Verbosity::get() > 1) std::cout << "ok\n";
	}

}
