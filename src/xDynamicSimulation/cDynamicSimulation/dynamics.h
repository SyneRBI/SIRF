/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <string>
#include <utility>
#include <vector>

#include "tissueparameters.h"
#include "tissuelabelmapper.h"


#include "gadgetron_data_containers.h"



typedef float TimeAxisType;
typedef float SignalAxisType;


typedef std::tuple<SignalAxisType,SignalAxisType,SignalAxisType> SignalBin;
typedef std::vector< std::pair<TimeAxisType, SignalAxisType> > SignalContainer;

// typedef ... MotionFieldContainer;

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

	std::vector<sirf::AcquisitionsVector> get_binned_mr_acquisitions( void )
	{
		return this->binned_mr_acquisitions_;
	};

	sirf::AcquisitionsVector get_binned_mr_acquisitions( int const bin_num)
	{
		if(bin_num >= this->num_simul_states_)
			throw std::runtime_error("Please access only bin numbers in the range of 0 and num_simul_states_-1.");
		
		return this->binned_mr_acquisitions_[bin_num];
	};

	SignalAxisType linear_interpolate_signal(TimeAxisType time_point);
	void bin_mr_acquisitions( sirf::AcquisitionsVector all_acquisitions );

protected:

	int num_simul_states_;
	virtual void set_bins( int const num_bins );

	std::vector< SignalBin > signal_bins_;
	SignalContainer dyn_signal_; 

	std::vector<sirf::AcquisitionsVector> binned_mr_acquisitions_;
};




class ContrastDynamic : public aDynamic {

public:
	ContrastDynamic():aDynamic(){};
	ContrastDynamic(int const num_simul_states); 


	TissueParameterList get_interpolated_tissue_params(SignalAxisType signal);

	void add_dynamic_label(LabelType l) { this->list_cont_var_labels_.push_back(l);};

	void set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1);

protected:

	virtual void set_bins(int const num_bins);

	std::vector< LabelType > list_cont_var_labels_;
	std::pair< TissueParameter, TissueParameter > tissue_parameter_extremes_;

};


class MotionDynamic : public aDynamic {

public:
	MotionDynamic();
	MotionDynamic(int const num_simul_states);

	~MotionDynamic();

	int get_which_motion_dynamic_am_i();
	int get_num_total_motion_dynamics();

	std::string get_temp_folder_name();

	bool make_temp_folder();
	bool delete_temp_folder();

protected:

	virtual void set_bins( int const num_bins );
	std::string setup_tmp_folder_name( void );

	// MotionFieldContainer displacement_field_;
	const static std::string temp_folder_path_;
	std::string temp_folder_name_ ;




	static int num_total_motion_dynamics_;
	int which_motion_dynamic_am_i_;

};


