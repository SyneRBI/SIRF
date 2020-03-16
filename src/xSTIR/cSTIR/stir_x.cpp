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

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"
#include "stir/is_null_ptr.h"
#include "stir/error.h"

#include "sirf/STIR/stir_x.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

#ifdef STIR_USE_LISTMODEDATA
    typedef ListModeData LMD;
    typedef ListRecord LMR;
#else
    typedef CListModeData LMD;
    typedef CListRecord LMR;
#endif

float ListmodeToSinograms::get_time_at_which_prompt_rate_exceeds_threshold(const float threshold) const
{
    if (input_filename.empty())
        throw std::runtime_error("ListmodeToSinograms::get_time_at_which_prompt_rate_exceeds_threshold: Filename missing");

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
	const int num_rings =
		lm_data_ptr->get_scanner_ptr()->get_num_rings();
	const int num_detectors_per_ring =
		lm_data_ptr->get_scanner_ptr()->get_num_detectors_per_ring();


	//*********** Finally, do the real work

	CPUTimer timer;
	timer.start();

	double time_of_last_stored_event = 0;
	long num_stored_events = 0;
	Array<2, float> data_fan_sums(IndexRange2D(num_rings, num_detectors_per_ring));
	fan_sums_sptr.reset(new std::vector<Array<2, float> >);

	// go to the beginning of the binary data
	lm_data_ptr->reset();

	// TODO have to use lm_data_ptr->get_proj_data_info_sptr() once STIR PR 108 is merged
	max_ring_diff_for_fansums = 60;
	if (*lm_data_ptr->get_scanner_ptr() != Scanner(Scanner::Siemens_mMR))
	{
		warning("This is not mMR data. Assuming all possible ring differences are in the listmode file");
		max_ring_diff_for_fansums = lm_data_ptr->get_scanner_ptr()->get_num_rings() - 1;
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
	assert(num_rings > 0);
	assert(data_fan_sums.get_min_index() == 0);
	num_detectors_per_ring = data_fan_sums[0].get_length();
	assert(num_detectors_per_ring > 0);
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

void
ListmodeToSinograms::estimate_randoms_()
{
	PETAcquisitionDataInFile acq_temp(template_proj_data_name.c_str());
	shared_ptr<ProjData> template_projdata_ptr = acq_temp.data();
	std::string filename = output_filename_prefix + "_randoms" + "_f1g1d0b0.hs";
	randoms_sptr = acq_temp.new_acquisition_data(); // filename);
	ProjData& proj_data = *randoms_sptr->data();

	const int num_rings =
		template_projdata_ptr->get_proj_data_info_ptr()->get_scanner_ptr()->
		get_num_rings();
	const int num_detectors_per_ring =
		template_projdata_ptr->get_proj_data_info_ptr()->get_scanner_ptr()->
		get_num_detectors_per_ring();
	DetectorEfficiencies& efficiencies = *det_eff_sptr;

	{
		const ProjDataInfoCylindricalNoArcCorr * const proj_data_info_ptr =
			dynamic_cast<const ProjDataInfoCylindricalNoArcCorr * const>
			(proj_data.get_proj_data_info_ptr());
		if (proj_data_info_ptr == 0)
		{
			error("Can only process not arc-corrected data\n");
		}

		const int mashing_factor =
			proj_data_info_ptr->get_view_mashing_factor();

		shared_ptr<Scanner>
			scanner_sptr(new Scanner(*proj_data_info_ptr->get_scanner_ptr()));
		unique_ptr<ProjDataInfo> uncompressed_proj_data_info_uptr
			(ProjDataInfo::construct_proj_data_info(scanner_sptr,
				/*span=*/1, max_ring_diff_for_fansums,
				/*num_views=*/num_detectors_per_ring / 2,
				scanner_sptr->get_max_num_non_arccorrected_bins(),
				/*arccorrection=*/false));
		const ProjDataInfoCylindricalNoArcCorr * const uncompressed_proj_data_info_ptr =
			dynamic_cast<const ProjDataInfoCylindricalNoArcCorr * const>
			(uncompressed_proj_data_info_uptr.get());
		Bin bin;
		Bin uncompressed_bin;

		for (bin.segment_num() = proj_data.get_min_segment_num();
			bin.segment_num() <= proj_data.get_max_segment_num();
			++bin.segment_num())
		{

			for (bin.axial_pos_num() = proj_data.get_min_axial_pos_num
				(bin.segment_num());
				bin.axial_pos_num() <= proj_data.get_max_axial_pos_num
				(bin.segment_num());
			++bin.axial_pos_num())
			{
				Sinogram<float> sinogram = proj_data_info_ptr->get_empty_sinogram
					(bin.axial_pos_num(), bin.segment_num());
				const float out_m = proj_data_info_ptr->get_m(bin);
				const int in_min_segment_num =
					proj_data_info_ptr->get_min_ring_difference(bin.segment_num());
				const int in_max_segment_num =
					proj_data_info_ptr->get_max_ring_difference(bin.segment_num());

				// now loop over uncompressed detector-pairs
				{
					for (uncompressed_bin.segment_num() = in_min_segment_num;
						uncompressed_bin.segment_num() <= in_max_segment_num;
						++uncompressed_bin.segment_num())
						for (uncompressed_bin.axial_pos_num() =
							uncompressed_proj_data_info_ptr->get_min_axial_pos_num
							(uncompressed_bin.segment_num());
					uncompressed_bin.axial_pos_num() <=
						uncompressed_proj_data_info_ptr->get_max_axial_pos_num
						(uncompressed_bin.segment_num());
					++uncompressed_bin.axial_pos_num())
						{
							const float in_m =
								uncompressed_proj_data_info_ptr->get_m(uncompressed_bin);
							if (fabs(out_m - in_m) > 1E-4)
								continue;

							// views etc
							if (proj_data.get_min_view_num() != 0)
								error("Can only handle min_view_num==0\n");
							for (bin.view_num() = proj_data.get_min_view_num();
								bin.view_num() <= proj_data.get_max_view_num();
								++bin.view_num())
							{

								for (bin.tangential_pos_num() = proj_data_info_ptr->get_min_tangential_pos_num();
									bin.tangential_pos_num() <= proj_data_info_ptr->get_max_tangential_pos_num();
									++bin.tangential_pos_num())
								{
									uncompressed_bin.tangential_pos_num() =
										bin.tangential_pos_num();
									for (uncompressed_bin.view_num() =
										bin.view_num()*mashing_factor;
										uncompressed_bin.view_num() <
										(bin.view_num() + 1)*mashing_factor;
									++uncompressed_bin.view_num())
									{
										int ra = 0, a = 0;
										int rb = 0, b = 0;
										uncompressed_proj_data_info_ptr->get_det_pair_for_bin(a, ra, b, rb,
											uncompressed_bin);
										/*(*segment_ptr)[bin.axial_pos_num()]*/
										sinogram[bin.view_num()][bin.tangential_pos_num()] +=
											efficiencies[ra][a] * efficiencies[rb][b%num_detectors_per_ring];
									}
								}
							}
						}
				}
				proj_data.set_sinogram(sinogram);
			}
		}
	}
	randoms_sptr->write(filename.c_str());
}

PETAcquisitionSensitivityModel::
PETAcquisitionSensitivityModel(PETAcquisitionData& ad)
{
	shared_ptr<PETAcquisitionData>
		sptr_ad(ad.new_acquisition_data());
	sptr_ad->inv(MIN_BIN_EFFICIENCY, ad);
	shared_ptr<BinNormalisation> 
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
	shared_ptr<BinNormalisation>
		sptr_n(new BinNormalisationFromECAT8(filename));
	//shared_ptr<BinNormalisation> sptr_0;
	//norm_.reset(new ChainedBinNormalisation(sptr_n, sptr_0));
	norm_ = sptr_n;
}

Succeeded 
PETAcquisitionSensitivityModel::set_up(const shared_ptr<ProjDataInfo>& sptr_pdi)
{
	return norm_->set_up(sptr_pdi);
}

void
PETAcquisitionSensitivityModel::unnormalise(PETAcquisitionData& ad) const
{
	BinNormalisation* norm = norm_.get();
	norm->undo(*ad.data(), 0, 1);
}

void
PETAcquisitionSensitivityModel::normalise(PETAcquisitionData& ad) const
{
	BinNormalisation* norm = norm_.get();
	norm->apply(*ad.data(), 0, 1);
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
PETAttenuationModel::unnormalise(PETAcquisitionData& ad) const
{
	//std::cout << "in PETAttenuationModel::unnormalise\n";
	BinNormalisation* norm = norm_.get();
	shared_ptr<DataSymmetriesForViewSegmentNumbers>
		symmetries_sptr(sptr_forw_projector_->get_symmetries_used()->clone());
	norm->undo(*ad.data(), 0, 1, symmetries_sptr);
}

void
PETAttenuationModel::normalise(PETAcquisitionData& ad) const
{
	BinNormalisation* norm = norm_.get();
	shared_ptr<DataSymmetriesForViewSegmentNumbers>
		symmetries_sptr(sptr_forw_projector_->get_symmetries_used()->clone());
	norm->apply(*ad.data(), 0, 1, symmetries_sptr);
}

//void
//PETAcquisitionModel::set_bin_efficiency
//(shared_ptr<PETAcquisitionData> sptr_data)
//{
//	shared_ptr<PETAcquisitionData>
//		sptr_ad(sptr_data->new_acquisition_data());
//	sptr_ad->inv(MIN_BIN_EFFICIENCY, *sptr_data);
//	sptr_normalisation_.reset
//		(new BinNormalisationFromProjData(sptr_ad->data()));
//	sptr_normalisation_->set_up(sptr_ad->get_proj_data_info_sptr());
//}

Succeeded 
PETAcquisitionModel::set_up(
	shared_ptr<PETAcquisitionData> sptr_acq,
	shared_ptr<STIRImageData> sptr_image)
{
	Succeeded s = Succeeded::no;
	if (sptr_projectors_.get()) {
		s = sptr_projectors_->set_up
			(sptr_acq->get_proj_data_info_sptr(), sptr_image->data_sptr());
		sptr_acq_template_ = sptr_acq;
		sptr_image_template_ = sptr_image;
	}
	if (s == Succeeded(Succeeded::yes)) {
		if (sptr_asm_ && sptr_asm_->data())
			s = sptr_asm_->set_up(sptr_acq->get_proj_data_info_sptr());
	}
	return s;
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
PETAcquisitionModel::forward(PETAcquisitionData& ad, const STIRImageData& image,
	int subset_num, int num_subsets, bool zero)
{
	shared_ptr<ProjData> sptr_fd = ad.data();
	sptr_projectors_->get_forward_projector_sptr()->forward_project
		(*sptr_fd, image.data(), subset_num, num_subsets, zero);

	float one = 1.0;

	if (sptr_add_.get()) {
		std::cout << "additive term added...";
		ad.axpby(&one, ad, &one, *sptr_add_);
		//ad.axpby(1.0, ad, 1.0, *sptr_add_);
		std::cout << "ok\n";
	}
	else
		std::cout << "no additive term added\n";

	PETAcquisitionSensitivityModel* sm = sptr_asm_.get();
	if (sm && sm->data() && !sm->data()->is_trivial()) {
		std::cout << "applying unnormalisation...";
		sptr_asm_->unnormalise(ad);
		std::cout << "ok\n";
	}
	else
		std::cout << "no unnormalisation applied\n";

	if (sptr_background_.get()) {
		std::cout << "background term added...";
		ad.axpby(&one, ad, &one, *sptr_background_);
		//ad.axpby(1.0, ad, 1.0, *sptr_background_);
		std::cout << "ok\n";
	}
	else
		std::cout << "no background term added\n";
}

shared_ptr<PETAcquisitionData>
PETAcquisitionModel::forward(const STIRImageData& image, 
	int subset_num, int num_subsets)
{
	shared_ptr<PETAcquisitionData> sptr_ad;
	sptr_ad = sptr_acq_template_->new_acquisition_data();
	shared_ptr<ProjData> sptr_fd = sptr_ad->data();
	//if (num_subsets > 1)
	//	sptr_fd->fill(0.0f);
	forward(*sptr_ad, image, subset_num, num_subsets, num_subsets > 1);

	//sptr_projectors_->get_forward_projector_sptr()->forward_project
	//	(*sptr_fd, image.data(), subset_num, num_subsets);
	////sptr_fd->fill(1.0f);

	//if (sptr_add_.get()) {
	//	std::cout << "additive term added...";
	//	sptr_ad->axpby(1.0, *sptr_ad, 1.0, *sptr_add_);
	//	std::cout << "ok\n";
	//}
	//else
	//	std::cout << "no additive term added\n";

	////if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
	//PETAcquisitionSensitivityModel* sm = sptr_asm_.get();
	//if (sm && sm->data() && !sm->data()->is_trivial()) {
	//	std::cout << "applying unnormalisation...";
	//	sptr_asm_->unnormalise(*sptr_ad);
	//	//sptr_normalisation_->undo(*sptr_fd, 0, 1);
	//	std::cout << "ok\n";
	//}
	//else
	//	std::cout << "no unnormalisation applied\n";

	//if (sptr_background_.get()) {
	//	std::cout << "background term added...";
	//	sptr_ad->axpby(1.0, *sptr_ad, 1.0, *sptr_background_);
	//	std::cout << "ok\n";
	//}
	//else
	//	std::cout << "no background term added\n";

	return sptr_ad;
}

shared_ptr<STIRImageData> 
PETAcquisitionModel::backward(PETAcquisitionData& ad, 
	int subset_num, int num_subsets)
{
	shared_ptr<STIRImageData> sptr_id;
	sptr_id = sptr_image_template_->new_image_data();
	shared_ptr<Image3DF> sptr_im = sptr_id->data_sptr();

	//if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
	PETAcquisitionSensitivityModel* sm = sptr_asm_.get();
	if (sm && sm->data() && !sm->data()->is_trivial()) {
		std::cout << "applying unnormalisation...";
		shared_ptr<PETAcquisitionData> sptr_ad(ad.new_acquisition_data());
		sptr_ad->fill(ad);
		sptr_asm_->unnormalise(*sptr_ad);
		//sptr_normalisation_->undo(*sptr_ad->data(), 0, 1);
		std::cout << "ok\n";
		std::cout << "backprojecting...";
		sptr_projectors_->get_back_projector_sptr()->back_project
			(*sptr_im, *sptr_ad, subset_num, num_subsets);
		std::cout << "ok\n";
	}
	else {
		std::cout << "backprojecting...";
		sptr_projectors_->get_back_projector_sptr()->back_project
			(*sptr_im, ad, subset_num, num_subsets);
		std::cout << "ok\n";
	}

	return sptr_id;
}
