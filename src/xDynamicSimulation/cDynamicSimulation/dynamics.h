/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include <utility>
#include <vector>


#include "gadgetron_data_containers.h"



typedef float TimeAxisType;
typedef float SignalAxisType;


typedef std::tuple<SignalAxisType,SignalAxisType,SignalAxisType> SignalBin;
typedef std::vector< std::pair<TimeAxisType, SignalAxisType> > SignalContainer;

// typedef ... MotionFieldContainer;

bool is_in_bin( SignalAxisType const signal, SignalBin const bin);


class aDynamic {

public:

	aDynamic() {};
	aDynamic(int const num_simul_states);


	std::vector< SignalBin > get_bins( void ){ return this->signal_bins_;};
	void set_num_simul_states(int const num_states);
	void set_dyn_signal(SignalContainer signal);

	SignalAxisType linear_interpolate_signal(TimeAxisType time_point);
	std::vector<AcquisitionsVector> bin_mr_acquisitions( AcquisitionsVector all_acquisitions );

protected:

	int num_simul_states_;
	void set_bins( int const num_bins );

	std::vector< SignalBin > signal_bins_;
	SignalContainer dyn_signal_; 


};




/*
class MotionDynamic : public aDynamic {

public:

protected:

	MotionFieldContainer motion_field_;

};
*/




