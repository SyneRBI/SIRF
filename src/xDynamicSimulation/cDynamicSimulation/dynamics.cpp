/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "dynamics.h"


aDynamic::aDynamic(int const num_simul_states) : num_simul_states_(num_simul_states)
{
	set_bins( num_simul_states );
}


void aDynamic::set_bins(int const num_bins)
{

	this->signal_bins_.clear();

	for(int i_state=0; i_state<num_bins; i_state++)
	{	
		SignalBin bin;

		std::get<0>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) - 1.f/(2*num_bins);
		std::get<1>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins);
		std::get<2>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) + 1.f/(2*num_bins);
		
		if( std::get<0>(bin) < 0 )
			std::get<0>(bin) = ( 1 + std::get<0>(bin) );

		if( std::get<1>(bin) < 0 )
			std::get<1>(bin) = ( 1 + std::get<1>(bin) );

		if( std::get<2>(bin) < 0 )
			std::get<2>(bin) = ( 1 + std::get<2>(bin) );

		this->signal_bins_.push_back( bin );
	}
}


void aDynamic::set_num_simul_states(int const num_states)
{
	this->num_simul_states_ = num_states;
	set_bins( num_states );

}

void aDynamic::set_dyn_signal(SignalContainer signal) 
{
	this->dyn_signal_ = signal;
}


SignalAxisType aDynamic::linear_interpolate_signal(TimeAxisType time_point)
{

	size_t const num_sig_points = this->dyn_signal_.size();

	size_t first_bigger_thant_time_point=-1;

	for( size_t i=0; i<num_sig_points; i++)
	{
		if(this->dyn_signal_[i].first > time_point)
		{
			first_bigger_thant_time_point = i;
			break;
		}
	}

	SignalAxisType interpol_signal;

	if( first_bigger_thant_time_point == 0)
		interpol_signal = this->dyn_signal_[0].second;
	else if( first_bigger_thant_time_point == -1 )
		interpol_signal = this->dyn_signal_[num_sig_points-1].second;
	else
	{
		interpol_signal = dyn_signal_[first_bigger_thant_time_point-1].second + 
							(time_point - dyn_signal_[first_bigger_thant_time_point-1].first )
						   *(dyn_signal_[first_bigger_thant_time_point].second - dyn_signal_[first_bigger_thant_time_point-1].second)
						   /(dyn_signal_[first_bigger_thant_time_point].first  - dyn_signal_[first_bigger_thant_time_point-1].first );
	}

	return interpol_signal;

}
