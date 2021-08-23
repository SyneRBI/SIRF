/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <string>
#include <utility>
#include <vector>
#include <stdexcept>

#include <ismrmrd/ismrmrd.h>

#include "sirf/cDynamicSimulation/phantom_input.h"
#include "sirf/cDynamicSimulation/tissueparameters.h"
#include "sirf/cDynamicSimulation/tissuelabelmapper.h"


#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "sirf/Reg/ImageWeightedMean.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"


#define SIRF_SCANNER_MS_PER_TIC 2.5


typedef float TimeAxisType;
typedef float SignalAxisType;


typedef std::tuple<SignalAxisType,SignalAxisType,SignalAxisType> SignalBin;
typedef std::pair<TimeAxisType, SignalAxisType> SignalPoint;
typedef std::vector< SignalPoint > SignalContainer;

// typedef std::vector< ISMRMRD::Image< DataTypeMotionFields > > MotionFieldContainer;
typedef sirf::NiftiImageData3DDisplacement<float> MotionFieldType;


bool is_in_bin( SignalAxisType const signal, SignalBin const bin);

typedef sirf::AcquisitionsVector MRDataType;
MRDataType intersect_mr_acquisition_data( const sirf::MRAcquisitionData& one_dat, const sirf::MRAcquisitionData& other_dat );



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

	SignalAxisType linear_interpolate_signal(const TimeAxisType time_point) const;

private:

	SignalContainer signal_; 
};

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

	void set_bins( void ){

		this->signal_bins_.clear();

		this->cyclic_ ? set_cyclic_bins(num_bins_) 
					  : set_non_cyclic_bins(num_bins_);
	}


private:
	void set_cyclic_bins( int const num_bins);
	void set_non_cyclic_bins( int const num_bins);

	int num_bins_;
	bool cyclic_;
	std::vector< SignalBin > signal_bins_;

};


class Dynamic{

public:
	
	Dynamic() : bp_(){}
	Dynamic(const unsigned int num_states) : bp_(num_states){
	}
	
	int get_num_simul_states(void) const {
		return bp_.get_num_bins();
	}

	void set_dynamic_signal(const SignalContainer& sig)
	{
		sp_.set_signal(sig);
	}

	std::vector< SignalBin > get_bins(void) const {
		 return bp_.get_bins();
	}

	void set_cyclicality(const bool cyclic){
		this->bp_.set_cyclicality(cyclic);
	}

	SignalAxisType interpolate_signal(const TimeAxisType time_point) const
	{
		return this->sp_.linear_interpolate_signal(time_point);
	}

protected: 
	SurrogateProcessor sp_;
	BinProcessor bp_;
};


class ContrastProcessor {

public:

	ContrastProcessor(){};
	
	TissueParameterList get_interpolated_tissue_params(SignalAxisType signal) const;
	
	void set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1);

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


class MotionProcessor {

public:
	MotionProcessor();
	~MotionProcessor();

	std::string get_temp_folder_name(void) const
	{
		return this->temp_folder_name_;
	}
	void set_ground_truth_folder_name(const std::string name_existing_folder_prefix)
	{
		this->ground_truth_folder_name_ = name_existing_folder_prefix;
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
	void prep_displacement_fields( void );
	// void align_motion_fields_with_image( const sirf::STIRImageData& img);

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

class MRDynamic : public Dynamic{

public:

	MRDynamic(unsigned int const num_simul_states): Dynamic(num_simul_states){}

	std::vector<sirf::AcquisitionsVector> get_binned_mr_acquisitions(void) const
	{
		return this->binned_mr_acquisitions_;
	}

	sirf::AcquisitionsVector get_binned_mr_acquisitions(unsigned int const bin_num) const
	{
		if(bin_num >= this->bp_.get_num_bins())
			throw std::runtime_error("Please access only bin numbers in the range of 0 and num_simul_states_-1.");
		
		return this->binned_mr_acquisitions_[bin_num];
	}

	virtual void bin_mr_acquisitions(sirf::MRAcquisitionData& all_acquisitions)=0;
	
protected:
	std::vector<sirf::AcquisitionsVector> binned_mr_acquisitions_;

};


class MRMotionDynamic : public MRDynamic{

public:

	MRMotionDynamic(unsigned int const num_simul_states): MRDynamic(num_simul_states){}
	virtual void bin_mr_acquisitions( sirf::MRAcquisitionData& all_acquisitions );

	void set_ground_truth_folder_name(const std::string name_existing_folder_prefix)
	{
		mp_.set_ground_truth_folder_name(name_existing_folder_prefix);
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
	
	void save_ground_truth_displacements( const std::vector< SignalAxisType >& gt_signal_points) const
	{
		mp_.save_ground_truth_displacements(gt_signal_points,bp_.is_cyclic());
	}

	void save_ground_truth_displacements(void) const
	{
		mp_.save_ground_truth_displacements(bp_.get_bin_centers(), bp_.is_cyclic());
	}

	bool delete_temp_folder() const
	{
		mp_.delete_temp_folder();
	}

private:
	MotionProcessor mp_;
};

class MRContrastDynamic: public MRDynamic {

public:
	MRContrastDynamic(unsigned int const num_simul_states): MRDynamic(num_simul_states){};
	virtual void bin_mr_acquisitions( sirf::MRAcquisitionData& all_acquisitions );

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
	static std::vector<sirf::AcquisitionsVector> binned_mr_acquisitions_;

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

	void set_ground_truth_folder_name(const std::string name_existing_folder_prefix)
	{
		mp_.set_ground_truth_folder_name(name_existing_folder_prefix);
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
