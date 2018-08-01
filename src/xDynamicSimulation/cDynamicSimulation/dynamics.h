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

typedef std::vector< std::pair<TimeAxisType, SignalAxisType> > SignalContainer;

// typedef ... MotionFieldContainer;




class aDynamic {

public:

	aDynamic() {};
	aDynamic(int num_simul_states) : num_simul_states_(num_simul_states){};

	void set_num_simul_states(int const num_states);
	void set_dyn_signal(SignalContainer signal);

	SignalAxisType linear_interpolate_signal(TimeAxisType time_point);
	std::vector<AcquisitionsVector> bin_mr_acquisitions( AcquisitionsVector acq_vec );

protected:

	int num_simul_states_;
	SignalContainer dyn_signal_; 
};




/*
class MotionDynamic : public aDynamic {

public:

protected:

	MotionFieldContainer motion_field_;

};
*/


