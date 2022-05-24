/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 - 2022 Physikalisch-Technische Bundesanstalt (PTB)

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
\ingroup Simulation
\brief Classes and utilities for handling dynamic processes in the simulation

\author Johannes Mayer
\author SyneRBI
*/


#pragma once

#include <string>
#include <utility>
#include <vector>
#include <stdexcept>

#include <ismrmrd/ismrmrd.h>

#include "sirf/cDynamicSimulation/phantom_input.h"
#include "sirf/cDynamicSimulation/tissueparameters.h"
#include "sirf/cDynamicSimulation/tissuelabelmapper.h"
#include "sirf/cDynamicSimulation/contrastgenerator.h"

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "sirf/Reg/ImageWeightedMean.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"


/*!
\brief TIC units used in ISMRMRD Encoding.Idx.TimeMr in miliseconds  
*/
#define SIRF_SCANNER_MS_PER_TIC 2.5

/*!
\brief Datatype used for the time axis of surrogate signals
*/
typedef float TimeAxisType;

/*!
\brief Datatype used for the signal axis of surrogate signals
*/
typedef float SignalAxisType;

/*!
\brief A signal bin consists of the triplet (lower border, centre, upper border)
*/
typedef std::tuple<SignalAxisType,SignalAxisType,SignalAxisType> SignalBin;

/*!
\brief Signal points of surrogate signals consist of a pair of (time,signal)
*/
typedef std::pair<TimeAxisType, SignalAxisType> SignalPoint; 

/*!
\brief Container class used to store surrogate signals
*/
typedef std::vector< SignalPoint > SignalContainer;

/*!
\brief A signal bin consists of the triplet (lower border, centre, upper border)
*/
typedef sirf::NiftiImageData3DDisplacement<float> MotionFieldType;

/*!
\brief Method to check if a signal point is inside a signal bin
*/
bool is_in_bin( SignalAxisType const signal, SignalBin const bin);

/*!
\brief The SIRF container class used to store MR acquisition data
*/
typedef sirf::AcquisitionsVector MRDataType;


/*!
\brief Method to intersect two MRAcquisitionData containers by compaing the acquisition scan_coutners.

Data are taken from one_dat, just scan_counters are compared. Acquisitions from other_dat are not returned!
*/
MRDataType intersect_mr_acquisition_data( const sirf::MRAcquisitionData& one_dat, const sirf::MRAcquisitionData& other_dat );


/*!
\brief Class to handle surrogate signals that are sets of points.
*/
class SurrogateProcessor{

public:
	SurrogateProcessor(){};

	void set_signal(const SignalContainer& sig){
		this->signal_ = sig;
	}

	SignalContainer get_signal(void) const{
		return signal_;
	}

	bool is_empty(void) const{
		return signal_.size() < 1;
	}

	SignalAxisType get_average_signal(const std::vector<TimeAxisType> timepoints) const{
		SignalAxisType avg_signal = 0.f;
		for(int i=0; i<timepoints.size(); ++i)
			avg_signal += linear_interpolate_signal(timepoints[i]);
		return avg_signal/timepoints.size();
	}
	/*!
	\brief Linear interpolate the signal from the discrete points held in signal_ any timepoint.
	Time points outside the time interval are extrapolated as constant. This assumes that the signal is sorted wrt. time points.
	*/
	SignalAxisType linear_interpolate_signal(const TimeAxisType time_point) const;

private:

	SignalContainer signal_; 
};


/*!
\brief Class to handle binning of surrogate signals.
*/
class BinProcessor{
public:

	BinProcessor(){
		num_bins_ = 1;
		cyclic_ = false;

		SignalBin bin;

		std::get<0>(bin) = 0.0;
		std::get<1>(bin) = 0.5;
		std::get<2>(bin) = 1.0;

		signal_bins_.push_back(bin);
	}

	BinProcessor(unsigned int const num_bins, bool const cyclic=false){
		
		if(num_bins < 1)
			throw std::runtime_error("Please a number of bins >= 1.");
		
		num_bins_ = num_bins;
		cyclic_ = cyclic;

		this->set_bins();
	}

	int get_num_bins( void ) const {
		return num_bins_;
	}
	
	std::vector<SignalBin> get_bins( void ) const {
		return signal_bins_;
	}

	std::vector<SignalAxisType> get_bin_centers(void) const {
		
		std::vector<SignalAxisType> centres;
		for(int i=0; i<signal_bins_.size(); ++i)
			centres.push_back(std::get<1>(signal_bins_.at(i)));
		
		return centres;
	}

	void set_cyclicality(const bool cyclic){
		this->cyclic_ = cyclic;
		set_bins();
	}
	bool is_cyclic(void) const{
		return cyclic_;
	}

	/*!
	\brief Overlays the interval [0,1] with num_bins_ non-overlapping bins
	*/
	void set_bins( void ){

		this->signal_bins_.clear();

		this->cyclic_ ? set_cyclic_bins(num_bins_) 
					  : set_non_cyclic_bins(num_bins_);
	}

private:
	/*!
	\brief Applies bins where the 0th bin is centered around 0 and assumes a cyclic interval [0,1].
	*/
	void set_cyclic_bins( int const num_bins);
	/*!
	\brief Applies bins where the 0th bin is centered around 0.5*1/num_bins_. 
	*/
	void set_non_cyclic_bins( int const num_bins);

	int num_bins_;
	bool cyclic_;
	std::vector< SignalBin > signal_bins_;

};

/*!
\brief Interface for simulated dynamic processes. 

Dynamics are generally composed of a SurrogateProcessor and a BinProcessor and a SurrogateProcessor and expose their methods.
*/
class Dynamic{

public:
	
	Dynamic() : bp_(){}
	Dynamic(const unsigned int num_states) : bp_(num_states){
	}
	
	int get_num_simul_states(void) const {
		return bp_.get_num_bins();
	}

	virtual void set_dynamic_signal(const SignalContainer& sig)
	{
		sp_.set_signal(sig);
	}

	std::vector< SignalBin > get_bins(void) const {
		 return bp_.get_bins();
	}

	virtual void set_cyclicality(const bool cyclic){
		this->bp_.set_cyclicality(cyclic);
	}

	SignalAxisType interpolate_signal(const TimeAxisType time_point) const
	{
		return this->sp_.linear_interpolate_signal(time_point);
	}

	/*!
	\brief Function to compute average surrogate signal over the timepoints at which MR rawdata was acquired

	The surrogate signal is interpolated linearly onto the timepoints at which the acquisitions were acquired.
	This subtracts the minimum time for the acquisition data assuming the first acquisition was acquires at t=0 ms.
	*/
	virtual SignalAxisType get_average_surrogate_signal(const sirf::MRAcquisitionData& ad) const
	{
		std::vector<TimeAxisType> timepts;
		for(int ia=0; ia<ad.number(); ++ia)
		{
			ISMRMRD::Acquisition acq;
			ad.get_acquisition(ia, acq);
			TimeAxisType acq_time_seconds = SIRF_SCANNER_MS_PER_TIC/1000.f * (TimeAxisType)acq.getHead().acquisition_time_stamp;
			timepts.push_back(acq_time_seconds);
		}
		TimeAxisType const t0 = *std::min_element(timepts.begin(), timepts.end());
		for(auto& element: timepts)
			element -= t0;
		
		return this->sp_.get_average_signal(timepts);
	}

protected: 
	SurrogateProcessor sp_;
	BinProcessor bp_;
};

/*!
\brief Class to handle the temporal change of tissue parameters to a surrogate signal.

Two tissue parameters can be passes as extreme points that correspond to a surrogate signal of 0 and 1.
E.g. two tissue parameters with different T1 due to inflow of a T1 contrast agent where 0 corresponds 
to the native T1 and 1 corresponds to the maximum observed concentration.

*/
class ContrastProcessor {

public:

	ContrastProcessor(){};
	
	TissueParameterList get_interpolated_tissue_params(SignalAxisType signal) const;
	
	void set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1);

	/*!
	/brief Add labels of a segmentation which behave identically (e.g. all labels containing blood)
	*/
	void add_dynamic_label(LabelType l) {
		this->list_cont_var_labels_.push_back(l);
	}

	std::vector< TimeAxisType > get_sampled_time_points(void) const {
		return time_points_sampled_;
	}

	void add_timepoint(const TimeAxisType time){
		time_points_sampled_.push_back(time);
	}

	void empty_timepoints(void){
		this->time_points_sampled_.clear();
	}

private:

	static std::vector< TimeAxisType > time_points_sampled_;

	std::vector< LabelType > list_cont_var_labels_;
	std::pair< TissueParameter, TissueParameter > tissue_parameter_extremes_;

};


/*!
\brief Utility class to process displacement fields.
*/
class MotionProcessor {

public:
	MotionProcessor();
	~MotionProcessor();

	std::string get_temp_folder_name(void) const
	{
		return this->temp_folder_name_;
	}
	void set_ground_truth_folder_name(const std::string prefix_output_dir)
	{
		this->ground_truth_folder_name_ = prefix_output_dir;
		this->make_ground_truth_folder();
	}

	void add_displacement_field(const MotionFieldType& dvf)
	{
		this->displacement_fields_.push_back( dvf );
	}

	int get_which_motion_processor_am_i(void) const
	{ 
		return this->which_motion_processor_am_i_; 
	}

	int get_num_total_motion_dynamics(void) const
	{ 
		return this->num_total_motion_processors_; 
	}

	sirf::NiftiImageData3DDeformation<float> get_interpolated_deformation_field(const SignalAxisType signal, const bool cyclic) const;

	void set_displacement_fields( const std::vector< MotionFieldType > &input_vectors);
	
	/*!
	/brief Function to store motion fields in a temporary folder instead of keeping them in memory.
	*/
	void prep_displacement_fields( void );
	

	void save_ground_truth_displacements( const std::vector< SignalAxisType >& gt_signal_points, const bool cyclic) const;

	bool delete_temp_folder() const;
	
protected:

	bool const destroy_upon_deletion_ = true;
	bool const keep_motion_fields_in_memory_ = true;

	std::string setup_tmp_folder_name( void );
	std::string setup_gt_folder_name( void );
	bool make_temp_folder() const;
	bool make_ground_truth_folder() const;

	sirf::NiftiImageData3DDeformation<float> calc_inverse_offset_deformation( sirf::NiftiImageData3DDeformation<float> offset_deformation ) const;
	sirf::NiftiImageData3DDisplacement<float> scale_displacementfields_to_mm( const sirf::NiftiImageData3DDisplacement<float> &dvf) const;

	std::string const temp_folder_prefix_  = "/media/sf_CCPPETMR/TestData/Input/xDynamicSimulation/cDynamicSimulation/Temp/";
	std::string const temp_mvf_prefix_ = "/motion_field_";

	std::string temp_folder_name_ ;
	std::string ground_truth_folder_name_ ;

	std::vector< std::string > temp_mvf_filenames_;

	static int num_total_motion_processors_;
	int which_motion_processor_am_i_;

	std::vector<MotionFieldType> displacement_fields_;
};

/*!
\brief Interface class to define dynamic processes for the simulation of MR acquisition. 
*/

class MRDynamic : public Dynamic{

public:
	MRDynamic(): Dynamic(){}
	MRDynamic(unsigned int const num_simul_states): Dynamic(num_simul_states){}

	std::vector<MRDataType> get_binned_mr_acquisitions(void) const
	{
		return this->binned_mr_acquisitions_;
	}

	MRDataType get_binned_mr_acquisitions(unsigned int const bin_num) const
	{
		if(bin_num >= binned_mr_acquisitions_.size())
			throw std::runtime_error("Please access only bin numbers in the range of 0 and num_simul_states_-1.");
		
		return binned_mr_acquisitions_.at(bin_num);
	}

	void clear_binning_data() 
	{
		std::vector<MRDataType> empty_acquis_vec;
		binned_mr_acquisitions_.swap( empty_acquis_vec );
		
		std::vector<std::deque<int> > empty_idx_corr;
		idx_corr_.swap(empty_idx_corr);
	
	}

	virtual std::deque<int> get_idx_corr(int const bin_num) const
	{
		if(bin_num >= binned_mr_acquisitions_.size())
			throw std::runtime_error("Please access only bin numbers in the range of 0 and num_simul_states_-1.");

		return this->idx_corr_.at(bin_num);
	}

	virtual std::vector<int> get_idx_corr_sizes() const
	{
		std::vector<int> idx_corr_sizes;
		for(int i=0; i<idx_corr_.size(); ++i)
			idx_corr_sizes.push_back(idx_corr_[i].size());

		return idx_corr_sizes;
	}

	virtual void bin_mr_acquisitions(sirf::MRAcquisitionData& all_acquisitions)=0;
	
protected:
	std::vector<MRDataType> binned_mr_acquisitions_;
	std::vector<std::deque<int> > idx_corr_;
};

/*!
\brief Implementation of MR dynamic executing motion during simulation. 
*/

class MRMotionDynamic : public MRDynamic{

public:

	MRMotionDynamic(unsigned int const num_simul_states): MRDynamic(num_simul_states){}
	virtual void bin_mr_acquisitions( sirf::MRAcquisitionData& all_acquisitions );

	void set_ground_truth_folder_name(const std::string prefix_output_dir)
	{
		mp_.set_ground_truth_folder_name(prefix_output_dir);
	}

	void add_displacement_field(const MotionFieldType& dvf)
	{
		mp_.add_displacement_field(dvf);
	}

	int get_which_motion_dynamic_am_i(void) const
	{ 
		return mp_.get_which_motion_processor_am_i();
	}

	int get_num_total_motion_dynamics(void) const
	{ 
		return mp_.get_num_total_motion_dynamics();
	}

	sirf::NiftiImageData3DDeformation<float> get_interpolated_deformation_field_at_timepoint(const TimeAxisType time_seconds) const
	{
		SignalAxisType sig = sp_.linear_interpolate_signal(time_seconds);			
		return this->get_interpolated_deformation_field(sig);
	}

	sirf::NiftiImageData3DDeformation<float> get_interpolated_deformation_field(const SignalAxisType signal) const
	{
		return mp_.get_interpolated_deformation_field(signal, bp_.is_cyclic());
	}

	/*!
	\brief Computes the average motion state over the time of acuqisition. 
	Subtracts the temporal offset i.e. assumes the first acquisition to be at t=0 ms
	*/
	sirf::NiftiImageData3DDeformation<float> get_average_deformation_field(const sirf::MRAcquisitionData& ad)
	{
		SignalAxisType const avg_sig = get_average_surrogate_signal(ad);
		return get_interpolated_deformation_field(avg_sig);
	}

	void set_displacement_fields( const std::vector< MotionFieldType > &input_vectors, bool const motion_fields_are_cyclic = false)
	{
		bp_.set_cyclicality(motion_fields_are_cyclic);
		mp_.set_displacement_fields(input_vectors);
	}
	void prep_displacement_fields(void)
	{
		mp_.prep_displacement_fields();
	}
	
	void save_ground_truth_displacements( const std::vector< SignalAxisType >& gt_signal_points) const
	{
		mp_.save_ground_truth_displacements(gt_signal_points, bp_.is_cyclic());
	}

	void save_ground_truth_displacements(void) const
	{
		this->save_ground_truth_displacements(bp_.get_bin_centers());
	}

	bool delete_temp_folder() const
	{
		mp_.delete_temp_folder();
	}


private:
	MotionProcessor mp_;
};

/*!
\brief Implementation of MR dynamic executing tissue parameter changes during simulation. 
*/

class MRContrastDynamic: public MRDynamic {

public:
	MRContrastDynamic(unsigned int const num_simul_states): MRDynamic(num_simul_states){};
	void bin_mr_acquisitions( sirf::MRAcquisitionData& all_acquisitions );

	TissueParameterList get_interpolated_tissue_params(SignalAxisType signal) const
	{
		return cp_.get_interpolated_tissue_params(signal);
	}
	
	void set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1)
	{
		cp_.set_parameter_extremes(tiss_at_0, tiss_at_1);
	}

	void add_dynamic_label(LabelType l) 
	{
		cp_.add_dynamic_label(l);
	}

	std::vector< TimeAxisType > get_sampled_time_points(void) const 
	{
		return cp_.get_sampled_time_points();
	}

private:
	ContrastProcessor cp_;
	
};


/*!
\brief Implementation of MR dynamic changing the signal according to a pre-computed magnetisation. 
*/

class ExternalMRContrastDynamic: public MRDynamic {

public:
	ExternalMRContrastDynamic(): MRDynamic(){};
	
	virtual void bin_mr_acquisitions( sirf::MRAcquisitionData& all_acquisitions )
	{
		if(external_signals_.size() != all_acquisitions.number())
		{
			std::stringstream message;
			message << "Please supply the same number (" << external_signals_.size() << ") of tissue signal sets in temporal order"
				"as number of readouts (" << all_acquisitions.number() << ") using set_tissue_signals() prior to calling this function.\n";
			throw std::runtime_error(message.str());
		}

		all_acquisitions.sort_by_time();

		vector<MRDataType>().swap(binned_mr_acquisitions_);

		ISMRMRD::Acquisition acq;
		for(int i=0; i<all_acquisitions.number(); ++i)
		{
			MRDataType av(all_acquisitions.acquisitions_info());	
			all_acquisitions.get_acquisition(i, acq);
			av.append_acquisition(acq);
			binned_mr_acquisitions_.push_back(av);
		}
	}

	void append_tissue_signals(const std::vector<ExternalTissueSignal>& ext_signals)
	{
		external_signals_.push_back(ext_signals);
	}

	void set_tissue_signals(std::vector<std::vector<ExternalTissueSignal> > signals)
	{
		external_signals_ = signals;
	}

	std::vector<ExternalTissueSignal> get_tissue_signals(const int i) const
	{
		return external_signals_[i];
	}
	
	int get_num_simul_states( void ) const
	{
		return binned_mr_acquisitions_.size();
	}

private:
	std::vector<std::vector<ExternalTissueSignal> > external_signals_;
};


// PET Dynamics and auxiliary methods

template <typename T>
struct Interval{

	Interval (): min_(0), max_(0) {}; 
	Interval (T min, T max): min_(min), max_(max) 
	{
		if( min > max )
			throw std::runtime_error("Please give an interval with smaller minimum than maximum.");
	};

	bool is_empty()
	{
		return std::abs(min_- max_) < 1e-8;
	}

	T min_;
	T max_;
};


template< typename T> 
Interval<T> intersect_intervals(const Interval<T>& one_interval, const Interval<T>& other_interval)
{
	Interval<T> intersection;

	if( one_interval.min_ > other_interval.max_ || other_interval.min_ > one_interval.max_ )
		return intersection;
	else 
	{
		intersection.min_ = std::max(one_interval.min_, other_interval.min_);
		intersection.max_ = std::min(one_interval.max_, other_interval.max_);

		return intersection;
	}
}

typedef Interval<TimeAxisType> TimeBin;
typedef std::vector< TimeBin > TimeBinSet;


TimeBin intersect_time_intervals( const TimeBin& one_interval, const TimeBin& other_interval);
TimeBinSet intersect_time_bin_sets( const TimeBinSet& one_set, const TimeBinSet& other_set);
TimeAxisType get_total_time_in_set( TimeBinSet& set_of_bins );
TimeAxisType get_time_from_between_two_signal_points(SignalAxisType signal, SignalPoint left_point, SignalPoint right_point);


class PETDynamic : public Dynamic{

public:

	PETDynamic(unsigned int const num_simul_states): Dynamic(num_simul_states){}

	void bin_total_time_interval(TimeBin time_interval_total_dynamic_process);

	TimeBinSet get_time_bin_set_for_state(unsigned int const which_state);
	TimeAxisType get_time_spent_in_bin(unsigned int const which_state );

protected:
	std::vector< TimeBinSet > binned_time_intervals_;

};


class PETMotionDynamic: public PETDynamic {

public:
	PETMotionDynamic(unsigned int const num_simul_states): PETDynamic(num_simul_states) {};

	void save_ground_truth_displacements(void) const
	{
		mp_.save_ground_truth_displacements(bp_.get_bin_centers(),bp_.is_cyclic());
	}
	
	void save_ground_truth_displacements( const std::vector< SignalAxisType >& gt_signal_points) const
	{
		mp_.save_ground_truth_displacements(gt_signal_points, bp_.is_cyclic());
	}

	void set_ground_truth_folder_name(const std::string prefix_output_dir)
	{
		mp_.set_ground_truth_folder_name(prefix_output_dir);
	}

	void add_displacement_field(const MotionFieldType& dvf)
	{
		mp_.add_displacement_field(dvf);
	}

	int get_which_motion_dynamic_am_i(void) const
	{ 
		return mp_.get_which_motion_processor_am_i();
	}

	int get_num_total_motion_dynamics(void) const
	{ 
		return mp_.get_num_total_motion_dynamics();
	}

	sirf::NiftiImageData3DDeformation<float> get_interpolated_deformation_field(const SignalAxisType signal) const
	{
		return mp_.get_interpolated_deformation_field(signal, bp_.is_cyclic());
	}

	void set_displacement_fields( const std::vector< MotionFieldType > &input_vectors, bool const motion_fields_are_cyclic = false)
	{
		bp_.set_cyclicality(motion_fields_are_cyclic);
		mp_.set_displacement_fields(input_vectors);
	}

	void prep_displacement_fields(void)
	{
		mp_.prep_displacement_fields();
	}
	
	bool delete_temp_folder() const
	{
		mp_.delete_temp_folder();
	}

private:
	MotionProcessor mp_;
	bool const keep_motion_fields_in_memory_ = true;
};

class PETContrastDynamic: public PETDynamic {

public:
	PETContrastDynamic(unsigned int const num_simul_states): PETDynamic(num_simul_states){};

	TissueParameterList get_interpolated_tissue_params(SignalAxisType signal) const
	{
		return cp_.get_interpolated_tissue_params(signal);
	}
	
	void set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1)
	{
		cp_.set_parameter_extremes(tiss_at_0, tiss_at_1);
	}

	void add_dynamic_label(LabelType l) 
	{
		cp_.add_dynamic_label(l);
	}

	std::vector< TimeAxisType > get_sampled_time_points(void) const 
	{
		return cp_.get_sampled_time_points();
	}

private:
	ContrastProcessor cp_;
};
