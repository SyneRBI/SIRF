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

#include "phantom_input.h"
#include "tissueparameters.h"
#include "tissuelabelmapper.h"


#include "gadgetron_data_containers.h"
#include "auxiliary_input_output.h"

#include "SIRFRegImageWeightedMean.h"
#include "SIRFImageDataDeformation.h"



typedef float TimeAxisType;
typedef float SignalAxisType;


typedef std::tuple<SignalAxisType,SignalAxisType,SignalAxisType> SignalBin;
typedef std::pair<TimeAxisType, SignalAxisType> SignalPoint;
typedef std::vector< SignalPoint > SignalContainer;

typedef std::vector< ISMRMRD::Image< DataTypeMotionFields > > MotionFieldContainer;



bool is_in_bin( SignalAxisType const signal, SignalBin const bin);

typedef sirf::AcquisitionsVector MRDataType;
MRDataType intersect_mr_acquisition_data(MRDataType one_dat, MRDataType other_dat);

class aDynamic {

public:

	aDynamic() {};
	aDynamic(int const num_simul_states);

	int get_num_simul_states( void ){ return this->num_simul_states_; };

	std::vector< SignalBin > get_bins( void ){ return this->signal_bins_;};
	void set_num_simul_states(int const num_states);
	void set_dyn_signal(SignalContainer signal);

	SignalAxisType linear_interpolate_signal(TimeAxisType time_point);
	
protected:

	int num_simul_states_;

	bool is_cyclic_dynamic_ = false;

	void set_bins( int const num_bins );
	void set_cyclic_bins( int const num_bins);
	void set_non_cyclic_bins( int const num_bins);


	std::vector< SignalBin > signal_bins_;
	SignalContainer dyn_signal_; 

};


class ContrastDynamic : virtual public aDynamic {

public:
	ContrastDynamic():aDynamic(){};
	ContrastDynamic(int const num_simul_states); 


	TissueParameterList get_interpolated_tissue_params(SignalAxisType signal);

	void add_dynamic_label(LabelType l) { this->list_cont_var_labels_.push_back(l);};

	void set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1);

protected:

	std::vector< LabelType > list_cont_var_labels_;
	std::pair< TissueParameter, TissueParameter > tissue_parameter_extremes_;

};


class MotionDynamic : virtual public aDynamic {

public:
	MotionDynamic();
	MotionDynamic(int const num_simul_states);

	~MotionDynamic();

	SIRFImageDataDeformation get_interpolated_displacement_field(SignalAxisType signal);

	int get_which_motion_dynamic_am_i();
	int get_num_total_motion_dynamics();

	std::string get_temp_folder_name();

	void set_displacement_fields( ISMRMRD::NDArray< DataTypeMotionFields >& motion_fields, bool const motion_fields_are_cyclic = false);
	     
	virtual void prep_displacements_fields( void );

	bool delete_temp_folder();

protected:

	bool const destroy_upon_deletion_ = true;
	bool const keep_motion_fields_in_memory_ = false;

	std::string setup_tmp_folder_name( void );
	bool make_temp_folder();

	std::string const temp_folder_prefix_  = "/tmp/";;
	std::string const temp_mvf_prefix_ = "/motion_field_";
	std::string temp_folder_name_ ;
	
	std::vector< std::string > temp_mvf_filenames_;

	static int num_total_motion_dynamics_;
	int which_motion_dynamic_am_i_;

	MotionFieldContainer displacment_fields_;
	std::vector< SIRFImageDataDeformation > sirf_displacement_fields_; 

};

class aMRDynamic : virtual public aDynamic{

public:

	aMRDynamic();
	aMRDynamic(int const num_simul_states);

	std::vector<sirf::AcquisitionsVector> get_binned_mr_acquisitions( void );
	sirf::AcquisitionsVector get_binned_mr_acquisitions( int const bin_num );
	void bin_mr_acquisitions( sirf::AcquisitionsVector all_acquisitions );

protected:

	std::vector<sirf::AcquisitionsVector> binned_mr_acquisitions_;


};



class MRMotionDynamic : public aMRDynamic, public MotionDynamic {


public:
	MRMotionDynamic():aMRDynamic(), MotionDynamic() {};
	MRMotionDynamic(int const num_simul_states): aMRDynamic(num_simul_states), MotionDynamic(num_simul_states) {};

	void prep_displacements_fields( void );
};

class MRContrastDynamic: public aMRDynamic, public ContrastDynamic {


public:
	MRContrastDynamic():aMRDynamic(), ContrastDynamic() {};
	MRContrastDynamic(int const num_simul_states): aMRDynamic(num_simul_states), ContrastDynamic(num_simul_states) {};
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


class aPETDynamic : virtual public aDynamic{

public:

	aPETDynamic();
	aPETDynamic(int const num_simul_states);

	void bin_total_time_interval(TimeBin time_interval_total_dynamic_process);

	TimeBinSet get_time_bin_set_for_state(unsigned int const which_state);
	TimeAxisType get_time_spent_in_bin(unsigned int const which_state );

protected:

	std::vector< TimeBinSet > binned_time_intervals_;

};


class PETMotionDynamic: public aPETDynamic, public MotionDynamic{

public:
	PETMotionDynamic():aPETDynamic(), MotionDynamic() {};
	PETMotionDynamic(int const num_simul_states): aPETDynamic(num_simul_states), MotionDynamic(num_simul_states) {};

	void align_motion_fields_with_image( const sirf::PETImageData& img);
	void prep_displacements_fields( void );
};

class PETContrastDynamic: public aPETDynamic, public ContrastDynamic {

public:
	PETContrastDynamic():aPETDynamic(), ContrastDynamic() {};
	PETContrastDynamic(int const num_simul_states): aPETDynamic(num_simul_states), ContrastDynamic(num_simul_states) {};
};
