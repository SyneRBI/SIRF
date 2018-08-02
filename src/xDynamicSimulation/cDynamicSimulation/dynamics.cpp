/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include <stdexcept>
#include <deque>

#include <ismrmrd/ismrmrd.h>


#include "dynamics.h"

bool is_in_bin( SignalAxisType const signal, SignalBin const bin)
{

	auto bin_min = std::get<0>(bin);
	auto bin_max = std::get<2>(bin);

	if( bin_min < bin_max )
		return (signal >= bin_min && signal < bin_max);
	else if ( bin_min > bin_max )
		return (signal >= bin_min || signal < bin_max);
	else
		return false;
}





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



void aDynamic::bin_mr_acquisitions( AcquisitionsVector all_acquisitions )
{

	if(this->dyn_signal_.size() == 0)
		throw std::runtime_error( "Please set a signal first. Otherwise you cannot bin your data, you dummy!" );


	std::deque< size_t > relevant_acq_numbers;

	for( size_t i=0; i<all_acquisitions.items(); i++)
		relevant_acq_numbers.push_back( i );

	for( int i_bin=0; i_bin<this->signal_bins_.size(); i_bin++)
	{

		auto bin = this->signal_bins_[i_bin];
		
		std::deque< size_t > acq_not_binned;

		AcquisitionsVector curr_acq_vector;
		curr_acq_vector.copy_acquisitions_info(all_acquisitions);
		
		ISMRMRD::Acquisition acq;

		while( relevant_acq_numbers.size() > 0 )	
		{
			acq_not_binned.clear();

			auto curr_pos = relevant_acq_numbers[0];
			std::cout << curr_pos << std::endl;

			relevant_acq_numbers.pop_front();	
			
			all_acquisitions.get_acquisition(curr_pos, acq);
			
			auto acq_hdr = acq.getHead();
			
			TimeAxisType acq_time = (TimeAxisType)acq_hdr.acquisition_time_stamp;
			std::cout << "nag " << std::endl;
			SignalAxisType signal_of_acq = this->linear_interpolate_signal( acq_time );
			std::cout << "nag " << std::endl;
			if( is_in_bin(signal_of_acq, bin) )
			{
				curr_acq_vector.append_acquisition(acq);
			}
			else
			{
				acq_not_binned.push_back(curr_pos);
			}
			std::cout << "nag " << std::endl;	
		}

		relevant_acq_numbers = acq_not_binned;
		this->binned_mr_acquisitions_.push_back( curr_acq_vector );
	}
}
























