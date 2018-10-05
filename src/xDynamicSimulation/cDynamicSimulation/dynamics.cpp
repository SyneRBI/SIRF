/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include <sstream>
#include <stdexcept>
#include <deque>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/system/error_code.hpp>

#include <ismrmrd/ismrmrd.h>


#include "dynamics.h"


using namespace sirf;


bool is_in_bin( SignalAxisType const signal, SignalBin const bin)
{

	auto bin_min = std::get<0>(bin);
	auto bin_max = std::get<2>(bin);

	if( bin_min < bin_max )
		return (signal >= bin_min && signal <= bin_max);
	else if ( bin_min > bin_max )
		return (signal >= bin_min || signal <= bin_max);
	else
		return false;
}


 
AcquisitionsVector intersect_mr_acquisition_data(AcquisitionsVector one_dat, AcquisitionsVector other_dat)
{

	bool one_dat_is_smaller = ( one_dat.items() >= other_dat.items() );

	typedef std::vector<uint32_t> CounterBox;

	CounterBox one_counters, other_counters;


	for( size_t i=0; i<one_dat.items(); i++)
	{
		ISMRMRD::Acquisition acq;

		one_dat.get_acquisition(i, acq);
		one_counters.push_back(acq.getHead().scan_counter);
	}

	for( size_t i=0; i<other_dat.items(); i++)
	{
		ISMRMRD::Acquisition acq;

		other_dat.get_acquisition(i, acq);
		other_counters.push_back(acq.getHead().scan_counter);
	}
	
	std::sort(one_counters.begin(), one_counters.end() );
	std::sort(other_counters.begin(), other_counters.end() );

	CounterBox intersected_counters(one_counters.size() + other_counters.size()); 

	CounterBox::iterator it;


	it = std::set_intersection( one_counters.begin(), one_counters.end(),
							    other_counters.begin(), other_counters.end(), 
							    intersected_counters.begin() );

	intersected_counters.resize( it - intersected_counters.begin() );
	

	MRDataType intersection;
	intersection.copy_acquisitions_info(one_dat);

	MRDataType& smaller_data_container = one_dat_is_smaller ? one_dat : other_dat;

	for( size_t i=0; i<smaller_data_container.items(); i++)
	{
		ISMRMRD::Acquisition acq;

		smaller_data_container.get_acquisition(i, acq);
		uint32_t acquis_counter = acq.getHead().scan_counter;
		if(std::find(intersected_counters.begin(), intersected_counters.end(), acquis_counter) != intersected_counters.end()) 
		{
			intersection.append_acquisition(acq);
    	} 
	}

	return intersection;

}


aDynamic::aDynamic(int const num_simul_states) : num_simul_states_(num_simul_states)
{
	set_bins( num_simul_states );
}

void aDynamic::set_bins(int const num_bins)
{
	this->signal_bins_.clear();

	this->signal_bins_.clear();


	if( this-> is_cyclic_dynamic_ )
		this->set_cyclic_bins(num_bins);
	else if ( ! this->is_cyclic_dynamic_)
		this->set_non_cyclic_bins(num_bins);
}

void aDynamic::set_cyclic_bins(int const num_bins)
{
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

void aDynamic::set_non_cyclic_bins(int const num_bins)
{
	for(int i_state=0; i_state<num_bins; i_state++)
	{	
		SignalBin bin;

		std::get<0>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins);
		std::get<1>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) + 1.f/(2*num_bins);
		std::get<2>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) + 1.f/(num_bins);
	
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





aMRDynamic::aMRDynamic(): aDynamic() {};
aMRDynamic::aMRDynamic(int const num_simul_states): aDynamic(num_simul_states){}

std::vector<sirf::AcquisitionsVector> aMRDynamic::get_binned_mr_acquisitions( void )
{
	return this->binned_mr_acquisitions_;
};

sirf::AcquisitionsVector aMRDynamic::get_binned_mr_acquisitions( int const bin_num )
{
	if(bin_num >= this->num_simul_states_)
		throw std::runtime_error("Please access only bin numbers in the range of 0 and num_simul_states_-1.");
	
	return this->binned_mr_acquisitions_[bin_num];
};

void aMRDynamic::bin_mr_acquisitions( AcquisitionsVector all_acquisitions )
{

	if(this->dyn_signal_.size() == 0)
		throw std::runtime_error( "Please set a signal first. Otherwise you cannot bin your data, you dummy!" );


	std::deque< size_t > relevant_acq_numbers;
	std::deque< size_t > acq_not_binned;


	for( size_t i=0; i<all_acquisitions.items(); i++)
		relevant_acq_numbers.push_back( i );

	for( int i_bin=0; i_bin<this->signal_bins_.size(); i_bin++)
	{

		auto bin = this->signal_bins_[i_bin];
	

		AcquisitionsVector curr_acq_vector;
		curr_acq_vector.copy_acquisitions_info(all_acquisitions);
		
		ISMRMRD::Acquisition acq;
		acq_not_binned.clear();

		while( relevant_acq_numbers.size() > 0 )	
		{
			auto curr_pos = relevant_acq_numbers[0];
			relevant_acq_numbers.pop_front();	
			
			all_acquisitions.get_acquisition(curr_pos, acq);
			
			auto acq_hdr = acq.getHead();
			
			TimeAxisType acq_time = (TimeAxisType)acq_hdr.acquisition_time_stamp;
			
			SignalAxisType signal_of_acq = this->linear_interpolate_signal( acq_time );
			
			if( is_in_bin(signal_of_acq, bin) )
			{
				curr_acq_vector.append_acquisition(acq);
			}
			else
			{
				acq_not_binned.push_back(curr_pos);
			}
			
		}
	
		relevant_acq_numbers.swap(acq_not_binned);
		this->binned_mr_acquisitions_.push_back( curr_acq_vector );
		
	}
}


ContrastDynamic::ContrastDynamic(int const num_simul_states) : aDynamic()
{ 
	this->num_simul_states_ =num_simul_states;
	this->set_bins(num_simul_states_);
}

void ContrastDynamic::set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1)
{
	this->tissue_parameter_extremes_.first = tiss_at_0;
	this->tissue_parameter_extremes_.second = tiss_at_1;
}


TissueParameterList ContrastDynamic::get_interpolated_tissue_params(SignalAxisType const signal)
{
	TissueParameterList tiss_list;

	for(size_t i=0; i< this->list_cont_var_labels_.size(); i++)
	{
		TissueParameter curr_par = ((1-signal) * this->tissue_parameter_extremes_.first + signal * this->tissue_parameter_extremes_.second);
		curr_par.name_ = "";	// name info is lost unfortunately, but better than the wrong information
		curr_par.label_ = this->list_cont_var_labels_[i];

		tiss_list.push_back( curr_par);
	}

	return tiss_list;
}


int MotionDynamic::num_total_motion_dynamics_ = 0;

MotionDynamic::MotionDynamic():aDynamic()
{
	this->which_motion_dynamic_am_i_ = num_total_motion_dynamics_;
	this->num_total_motion_dynamics_ += 1;
	
	this->temp_folder_name_ = setup_tmp_folder_name();

}

MotionDynamic::MotionDynamic(int const num_simul_states) : aDynamic()
{
	this->num_simul_states_ =num_simul_states;
	this->set_bins(num_simul_states_);

	this->which_motion_dynamic_am_i_ = num_total_motion_dynamics_;
	this->num_total_motion_dynamics_ += 1;

	this->temp_folder_name_ = setup_tmp_folder_name();

}


MotionDynamic::~MotionDynamic()
{ 
	// if( this->destroy_upon_deletion_)
	// 	this->delete_temp_folder();

	this->num_total_motion_dynamics_ -= 1; 
}


SIRFImageDataDeformation MotionDynamic::get_interpolated_displacement_field(SignalAxisType signal)
{
	if (signal > 1.f || signal< 0.f)
		throw std::runtime_error("Please pass a signal in the range of [0,1].");

	if( this->temp_displacement_fields_.size() == 0)
		throw std::runtime_error("Please use prep_displacements_fields() before calling get_interpolated_displacement_field");
	

	// check in which interval the signal lies
	SignalAxisType signal_on_bin_range;
	
	size_t const num_motion_fields = this->temp_mvf_filenames_.size();

	if( this->is_cyclic_dynamic_ )
		signal_on_bin_range = num_motion_fields * signal;
	else
		signal_on_bin_range = (num_motion_fields  - 1)* signal;

	int const bin_floor = int( signal_on_bin_range +1) -1;
	int const bin_ceil  = int( signal_on_bin_range + 1) % this->num_simul_states_;
	

	SignalAxisType const linear_interpolation_weight = signal_on_bin_range - bin_floor;

	  /// Constructor
    SIRFRegImageWeightedMean4D dvf_interpolator;
    
    dvf_interpolator.add_image( this->temp_displacement_fields_[bin_floor], 1 - linear_interpolation_weight);
    dvf_interpolator.add_image( this->temp_displacement_fields_[bin_ceil], linear_interpolation_weight);

    dvf_interpolator.update();
    
    SIRFImageDataDeformation output_deformation = dvf_interpolator.get_output();

    return output_deformation;

}




int MotionDynamic::get_which_motion_dynamic_am_i(){ return this->which_motion_dynamic_am_i_; }
int MotionDynamic::get_num_total_motion_dynamics(){ return this->num_total_motion_dynamics_; }

std::string MotionDynamic::setup_tmp_folder_name()
{
	std::string const current_folder_prefix = "temp_folder_motion_dyn_";
	std::stringstream tmp_stream;
	tmp_stream << this->temp_folder_prefix_ << current_folder_prefix << this->which_motion_dynamic_am_i_;
	return tmp_stream.str();

}

std::string MotionDynamic::get_temp_folder_name()
{
	return this->temp_folder_name_;
}

bool MotionDynamic::make_temp_folder()
{
	try
	{
		std::cout << "Generating temporary folder " << this->temp_folder_name_ << std::endl;
		boost::filesystem::path dir_to_make(this->temp_folder_name_.c_str());
		bool folder_creation_worked = boost::filesystem::create_directories(dir_to_make);
		return folder_creation_worked;

	}
	catch(boost::system::error_code& e)
	{
		std::cout << e.message() << std::endl;
		throw e;	
	}
}

bool MotionDynamic::delete_temp_folder()
{
	try
	{
		boost::filesystem::path dir_to_del( this->temp_folder_name_.c_str() );
		
		if( boost::filesystem::exists(dir_to_del) )
		{
			std::cout << "Deleting temporary folder " << this->temp_folder_name_ << std::endl;

			bool folder_deletion_worked = boost::filesystem::remove_all(dir_to_del);
			return folder_deletion_worked;
		}
		else
		{
			std::cout << "Folder " << this->temp_folder_name_ << " does not exist. Deletion omitted" << std::endl;
			return false;
		}

	}
	catch(boost::system::error_code& e)
	{
		std::cout << e.message() << std::endl;
		throw e;	
	}
;
}


void MotionDynamic::set_displacement_fields( ISMRMRD::NDArray< DataTypeMotionFields >& motion_fields, bool const motion_fields_are_cyclic)
{
	
	if ( motion_fields_are_cyclic )
	{
		this->is_cyclic_dynamic_ = true;
		this->set_bins( this->num_simul_states_ );
	}	


	using namespace ISMRMRD;

	const size_t* dimensions = motion_fields.getDims();

	size_t const Nt = dimensions[0];
	size_t const Nv = dimensions[1];
	size_t const Nz = dimensions[2];
	size_t const Ny = dimensions[3];
	size_t const Nx = dimensions[4];

	for(size_t nt=0; nt<Nt; nt++)
	{
		
		Image<DataTypeMotionFields> img(dimensions[4],dimensions[3], dimensions[2], dimensions[1]);
 		
 		for(uint16_t  nv= 0; nv<Nv ; nv++)
		for(uint16_t  nz= 0; nz<Nz ; nz++)
		for(uint16_t  ny= 0; ny<Ny ; ny++)
		for(uint16_t  nx= 0; nx<Nx ; nx++)
		{
			// size_t const lin_index = ((((Nt-1 -nt)*Nv + Nv-1 -nv)*Nz + Nz-1 - nz)*Ny + Ny-1 - ny)*Nx + Nx-1 - nx;
			size_t const lin_index = (((nt*Nv + nv)*Nz + nz)*Ny + ny)*Nx + nx;
			img(nx,ny,nz,nv) = 	  *(motion_fields.begin() + lin_index);
		}
		this->displacment_fields_.push_back(img);
	}
}

void MotionDynamic::prep_displacements_fields()
{
	if(this->displacment_fields_.size() == 0)
		throw std::runtime_error("Please call set_displacements_fields() first.");

	bool const temp_folder_creation_successful = this->make_temp_folder();

	if( temp_folder_creation_successful )
	{
		for(int i=0; i<this->displacment_fields_.size(); i++)
		{
			std::stringstream temp_filename_mvf;
			temp_filename_mvf << this->get_temp_folder_name() << this->temp_mvf_prefix_ << i;

			data_io::write_ISMRMRD_Image_to_Analyze<DataTypeMotionFields> (temp_filename_mvf.str(), this->displacment_fields_[i]);
			temp_filename_mvf << ".hdr";
			this->temp_mvf_filenames_.push_back(temp_filename_mvf.str());
		}

		if( this-> keep_motion_fields_in_memory_ == false)
			this->displacment_fields_ = MotionFieldContainer();
	}
	else
		throw std::runtime_error("The parent directory generation failed. Give a path to which thou hast access rights. Or maybe the directory already exists. This is dangerous. Then you should definitely choose a different temporary folder name.");

	for( size_t i=0; i<temp_mvf_filenames_.size(); i++)
	{
		SIRFImageDataDeformation temp_deformation( this->temp_mvf_filenames_[i] );
		temp_displacement_fields_.push_back( temp_deformation );
	}

	this->delete_temp_folder();
}




TimeBin intersect_time_intervals( const TimeBin& one_interval, const TimeBin& other_interval)
{
	return intersect_intervals<TimeAxisType>(one_interval, other_interval);
}

SetTimeBins intersect_set_time_bins( const SetTimeBins& one_set, const SetTimeBins& other_set)
{
	SetTimeBins intersected_set;
	for(size_t i=0; i<one_set.size();i++ )
	for(size_t j=0; j<other_set.size();j++ )
	{
		TimeBin temp_intersect = intersect_time_intervals(one_set[i], other_set[j]);
		if( !temp_intersect.is_empty() )			
			intersected_set.push_back( temp_intersect );
	}
	return intersected_set;
}

aPETDynamic::aPETDynamic(int const num_simul_states): aDynamic(num_simul_states){}


void aPETDynamic::bin_total_time_interval(TimeBin time_interval_total_dynamic_process)
{
	if(this->dyn_signal_.size() == 0)
		throw std::runtime_error( "Please set a signal first. Otherwise you cannot bin your data, you dummy!" );
	
	size_t const num_bins = signal_bins_.size();
	size_t const num_signal_supports = dyn_signal_.size();
	
	TimeAxisType leftmost_left_edge = std::min<TimeAxisType>( time_interval_total_dynamic_process.min_, dyn_signal_[0].first );
	TimeAxisType rightmost_right_edge = std::max<TimeAxisType>( time_interval_total_dynamic_process.max_,  dyn_signal_[num_signal_supports-1].first );

	for( size_t i_bin=0; i_bin<num_bins; i_bin++)
	{
		std::cout << "BIN #: " << i_bin <<std::endl;
		SignalBin bin = this->signal_bins_[i_bin];

		auto bin_min = std::get<0>(bin);
		auto bin_max = std::get<2>(bin);

		std::cout << "bin_min= "<< bin_min <<std::endl;
		std::cout << "bin_max= "<< bin_max <<std::endl;
	
		SetTimeBins time_intervals_for_bin;

		std::vector< TimeAxisType > left_bin_edges, right_bin_edges;

		for(size_t i_sig_pt=0; i_sig_pt<num_signal_supports-1; i_sig_pt++)
		{	
			SignalAxisType signal_this = dyn_signal_[i_sig_pt].second;
			SignalAxisType signal_next = dyn_signal_[i_sig_pt+1].second;
			
			bool const min_is_between_points = (bin_min >= signal_this && bin_min <= signal_next) || (bin_min <= signal_this && bin_min >= signal_next);
			bool const max_is_between_points = (bin_max >= signal_this && bin_max <= signal_next) || (bin_max <= signal_this && bin_max >= signal_next);
			
			TimeAxisType intersection_point;

			if( min_is_between_points )
			{
				intersection_point = get_time_from_between_two_signal_points(bin_min, dyn_signal_[i_sig_pt], dyn_signal_[i_sig_pt+1]);
				SignalAxisType f0 = this->linear_interpolate_signal(intersection_point);
				
				TimeAxisType delta_time = 0.1;
				SignalAxisType f_plus = this->linear_interpolate_signal( intersection_point + delta_time);
				SignalAxisType f_minus = this->linear_interpolate_signal(intersection_point - delta_time);

				if( f_plus > f0 && f_minus < f0)
					left_bin_edges.push_back(intersection_point);
				else if( f_plus < f0 && f_minus > f0)
					right_bin_edges.push_back(intersection_point);	
				

			}
					
			if( max_is_between_points )
			{
				intersection_point = get_time_from_between_two_signal_points(bin_max, dyn_signal_[i_sig_pt], dyn_signal_[i_sig_pt+1]);
				SignalAxisType f0 = this->linear_interpolate_signal(intersection_point);
				
				TimeAxisType delta_time = 0.1;
				SignalAxisType f_plus = this->linear_interpolate_signal( intersection_point + delta_time);
				SignalAxisType f_minus = this->linear_interpolate_signal(intersection_point - delta_time);

				if( f_plus > f0 && f_minus < f0)
					right_bin_edges.push_back(intersection_point);
				else if( f_plus < f0 && f_minus > f0)
					left_bin_edges.push_back(intersection_point);	
			}
		}
		std::sort( left_bin_edges.begin(), left_bin_edges.end() );
		std::sort( right_bin_edges.begin(), right_bin_edges.end() );

		for(size_t i=0; i<left_bin_edges.size(); i++)
			std::cout << "LBE " << i << " " << left_bin_edges[i] <<std::endl;

		for(size_t i=0; i<right_bin_edges.size(); i++)
			std::cout << "RBE " << i << " " << right_bin_edges[i] <<std::endl;

		size_t const num_left_edges = left_bin_edges.size();
		size_t const num_right_edges = right_bin_edges.size();

		if ( std::abs( num_left_edges - num_right_edges ) > 1 )
			throw std::runtime_error( "You got yourself a very weird combination of signal and bin number in your simulation. Consider passing another dynamic signal please.");	


		if( num_left_edges == num_right_edges )
		{
			if(left_bin_edges[0] < right_bin_edges[0])
			{
				for(size_t i=0; i<num_left_edges; i++)
				{
					TimeBin curr_bin(left_bin_edges[i], right_bin_edges[i]);
					time_intervals_for_bin.push_back(curr_bin);
				}
			}
			else
			{
				for(size_t i=0; i<num_left_edges-1; i++)
				{
					TimeBin curr_bin(left_bin_edges[i], right_bin_edges[i+1]);
					time_intervals_for_bin.push_back(curr_bin);
				}

				TimeBin leftmost_bin( leftmost_left_edge , right_bin_edges[0]);
				TimeBin rightmost_bin(left_bin_edges[num_left_edges-1], rightmost_right_edge);
				
				time_intervals_for_bin.push_back(leftmost_bin);
				time_intervals_for_bin.push_back(rightmost_bin);

			}
		}
		else if(num_left_edges > num_right_edges)	
		{
			for(size_t i=0; i<num_left_edges-1; i++)
			{
				TimeBin curr_bin(left_bin_edges[i], right_bin_edges[i]);
				time_intervals_for_bin.push_back(curr_bin);
			}
			TimeBin rightmost_bin( left_bin_edges[num_left_edges-1], rightmost_right_edge );
			time_intervals_for_bin.push_back(rightmost_bin);
		}
		else if(num_right_edges > num_left_edges)
		{
			for(size_t i=0; i<num_right_edges-1; i++)
			{
				TimeBin curr_bin(left_bin_edges[i], right_bin_edges[i+1]);
				time_intervals_for_bin.push_back(curr_bin);
			}
			TimeBin leftmost_bin( right_bin_edges[0], leftmost_left_edge );
			time_intervals_for_bin.push_back(leftmost_bin);

		}

		binned_time_intervals_.push_back(time_intervals_for_bin);
	}
}


TimeAxisType get_time_from_between_two_signal_points(SignalAxisType signal, SignalPoint left_point, SignalPoint right_point)
{
	if(std::abs(right_point.second - left_point.second) < 1e-8)
		return (right_point.first + left_point.first)/TimeAxisType(2);
	else
		return (signal-left_point.second) * (right_point.first - left_point.first) / (right_point.second - left_point.second) + left_point.first;
}




























