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


#include "sirf/cGadgetron/gadgetron_data_containers.h"
#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "sirf/cReg/ImageWeightedMean.h"
#include "sirf/cReg/NiftiImageData3DDeformation.h"


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
MRDataType intersect_mr_acquisition_data( MRDataType& one_dat, MRDataType& other_dat );

class aDynamic {

public:

	aDynamic() {};
	aDynamic(int const num_simul_states);

	virtual int get_num_simul_states( void ){ return this->num_simul_states_; };

	std::vector< SignalBin > get_bins( void ){ return this->signal_bins_;};
	void set_num_simul_states(int const num_states);
	void set_dyn_signal(SignalContainer signal);

	SignalAxisType linear_interpolate_signal(TimeAxisType time_point);
	
protected:

	int num_simul_states_;

	bool is_cyclic_dynamic_ = false;

	virtual void set_bins( int const num_bins );
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

	std::vector< TimeAxisType > get_time_points_sampled (void)const {return this->time_points_sampled_;};

	virtual int get_num_simul_states( void ){ return this->num_simul_states_; };


protected:

	static int num_simul_states_;
	static std::vector< TimeAxisType > time_points_sampled_;

	std::vector< LabelType > list_cont_var_labels_;
	std::pair< TissueParameter, TissueParameter > tissue_parameter_extremes_;
	virtual void set_bins( int const num_bins );

};


class MotionDynamic : virtual public aDynamic {

public:
	MotionDynamic();
	MotionDynamic(int const num_simul_states);

	~MotionDynamic();

	sirf::NiftiImageData3DDeformation<float> get_interpolated_deformation_field(SignalAxisType signal);

	int get_which_motion_dynamic_am_i();
	int get_num_total_motion_dynamics();

	std::string get_temp_folder_name();

	void set_ground_truth_folder_name( std::string const name_existing_folder );

	void set_displacement_fields( ISMRMRD::NDArray< DataTypeMotionFields >& motion_fields, bool const motion_fields_are_cyclic = false);
	void set_displacement_fields( std::vector< MotionFieldType > &input_vectors, bool const motion_fields_are_cyclic = false);
	     
	virtual void prep_displacement_fields( void );

	virtual void save_ground_truth_displacements( void );
	void save_ground_truth_displacements( std::vector< SignalAxisType > gt_signal_points);

	bool delete_temp_folder();
	
protected:

	bool const destroy_upon_deletion_ = true;
	bool const keep_motion_fields_in_memory_ = true;

	std::string setup_tmp_folder_name( void );
	std::string setup_gt_folder_name( void );
	bool make_temp_folder();
	bool make_ground_truth_folder();

	sirf::NiftiImageData3DDeformation<float> calc_inverse_offset_deformation( sirf::NiftiImageData3DDeformation<float> offset_deformation );

	sirf::NiftiImageData3DDisplacement<float> scale_displacementfields_to_mm( const sirf::NiftiImageData3DDisplacement<float> &dvf);

	std::string const temp_folder_prefix_  = "/media/sf_SharedFolder/CCPPETMR/PublicationData/Output/MotionFields/";
	std::string const temp_mvf_prefix_ = "/motion_field_";

	std::string temp_folder_name_ ;
	std::string ground_truth_folder_name_ ;

	std::vector< std::string > temp_mvf_filenames_;

	static int num_total_motion_dynamics_;
	int which_motion_dynamic_am_i_;

	std::vector<MotionFieldType> displacement_fields_;
	// std::vector< sirf::NiftiImageData3DDeformation<float> > sirf_displacement_fields_;
};

class aMRDynamic : virtual public aDynamic{

public:

	aMRDynamic();
	aMRDynamic(int const num_simul_states);

	virtual std::vector<sirf::AcquisitionsVector> get_binned_mr_acquisitions( void );
	virtual sirf::AcquisitionsVector get_binned_mr_acquisitions( int const bin_num );

	virtual void bin_mr_acquisitions( sirf::AcquisitionsVector& all_acquisitions )=0;

protected:

	std::vector<sirf::AcquisitionsVector> binned_mr_acquisitions_;


};



class MRMotionDynamic : public aMRDynamic, public MotionDynamic {


public:
	MRMotionDynamic():aMRDynamic(), MotionDynamic() {};
	MRMotionDynamic(int const num_simul_states): aMRDynamic(num_simul_states), MotionDynamic(num_simul_states) {};

	// void prep_displacement_fields( void );
	virtual void bin_mr_acquisitions( sirf::AcquisitionsVector& all_acquisitions );
};

class MRContrastDynamic: public aMRDynamic, public ContrastDynamic {


public:
	MRContrastDynamic():aMRDynamic(), ContrastDynamic() {};
	MRContrastDynamic(int const num_simul_states): aMRDynamic(num_simul_states), ContrastDynamic(num_simul_states) {};
	virtual void bin_mr_acquisitions( sirf::AcquisitionsVector& all_acquisitions );

	virtual std::vector<sirf::AcquisitionsVector> get_binned_mr_acquisitions( void );
	virtual sirf::AcquisitionsVector get_binned_mr_acquisitions( int const bin_num );


protected:
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

	void align_motion_fields_with_image( const sirf::STIRImageData& img);
	// void prep_displacement_fields( void );
private:
	bool const keep_motion_fields_in_memory_ = true;

};

class PETContrastDynamic: public aPETDynamic, public ContrastDynamic {

public:
	PETContrastDynamic():aPETDynamic(), ContrastDynamic() {};
	PETContrastDynamic(int const num_simul_states): aPETDynamic(num_simul_states), ContrastDynamic(num_simul_states) {};
};
