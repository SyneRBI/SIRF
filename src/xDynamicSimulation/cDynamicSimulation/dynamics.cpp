/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <deque>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/system/error_code.hpp>

#include <ismrmrd/ismrmrd.h>


#include "sirf/common/multisort.h"
#include "sirf/cDynamicSimulation/dynamics.h"

#include "sirf/Reg/NiftyResample.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"

#include <_reg_localTrans.h>

using namespace std;
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

// data are taken from one_dat, just scan_counters are compared
AcquisitionsVector intersect_mr_acquisition_data( const MRAcquisitionData& one_dat, const MRAcquisitionData& other_dat)
{
	typedef std::vector<uint32_t> CounterBox;

	CounterBox one_counters, other_counters;

	ISMRMRD::Acquisition acq;

	std::map<uint32_t, int> map_counter_idx;
	
	for( size_t i=0; i<one_dat.items(); i++)
	{
		one_dat.get_acquisition( i, acq );
		one_counters.push_back(acq.getHead().scan_counter);
		map_counter_idx[acq.getHead().scan_counter] = i;
	}

	for( size_t i=0; i<other_dat.items(); i++)
	{
		other_dat.get_acquisition( i, acq );
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

	for( int i=0; i<intersected_counters.size(); ++i)
	{
		int const idx_for_acq = map_counter_idx[intersected_counters[i]];
		one_dat.get_acquisition(idx_for_acq, acq);
		intersection.append_acquisition(acq);
	}

	return intersection;
}

void BinProcessor::set_cyclic_bins(int const num_bins)
{
	for(int i_state=0; i_state<num_bins; i_state++)
	{	
		SignalBin bin;

		std::get<0>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) - 1.f/(2*num_bins);
		std::get<1>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins);
		std::get<2>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) + 1.f/(2*num_bins);
		
		std::get<0>(bin) = std::get<0>(bin) < 0 ? ( 1 + std::get<0>(bin) ) : std::get<0>(bin);
		std::get<1>(bin) = std::get<1>(bin) < 0 ? ( 1 + std::get<1>(bin) ) : std::get<1>(bin);
		std::get<2>(bin) = std::get<2>(bin) < 0 ? ( 1 + std::get<2>(bin) ) : std::get<2>(bin);

		this->signal_bins_.push_back( bin );
	}
}

void BinProcessor::set_non_cyclic_bins(int const num_bins)
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


SignalAxisType SurrogateProcessor::linear_interpolate_signal(const TimeAxisType time_point) const
{

	size_t const num_sig_points = this->signal_.size();
	
	size_t first_bigger_thant_time_point=-1;

	for( size_t i=0; i<num_sig_points; i++)
	{
		if(this->signal_[i].first > time_point)
		{
			first_bigger_thant_time_point = i;
			break;
		}
	}

	SignalAxisType interpol_signal;

	if( first_bigger_thant_time_point == 0)
		interpol_signal = this->signal_[0].second;
	else if( first_bigger_thant_time_point == -1 )
		interpol_signal = this->signal_[num_sig_points-1].second;
	else
	{
		interpol_signal = signal_[first_bigger_thant_time_point-1].second + 
							(time_point - signal_[first_bigger_thant_time_point-1].first )
						   *(signal_[first_bigger_thant_time_point].second - signal_[first_bigger_thant_time_point-1].second)
						   /(signal_[first_bigger_thant_time_point].first  - signal_[first_bigger_thant_time_point-1].first );
	}

	return interpol_signal;

}

// static member variable initialization
std::vector< TimeAxisType > ContrastProcessor::time_points_sampled_ = std::vector<TimeAxisType>{};

void ContrastProcessor::set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1)
{
	this->tissue_parameter_extremes_.first = tiss_at_0;
	this->tissue_parameter_extremes_.second = tiss_at_1;
}


TissueParameterList ContrastProcessor::get_interpolated_tissue_params(SignalAxisType const signal) const
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


int MotionProcessor::num_total_motion_processors_ = 0;

MotionProcessor::MotionProcessor() 
{
	this->which_motion_processor_am_i_ = num_total_motion_processors_;
	this->num_total_motion_processors_ += 1;
	
	this->temp_folder_name_ = setup_tmp_folder_name();
	this->ground_truth_folder_name_ = setup_gt_folder_name();
}

MotionProcessor::~MotionProcessor()
{ 
	// if( this->destroy_upon_deletion_)
	// this->delete_temp_folder();

	this->num_total_motion_processors_ -= 1; 
}


NiftiImageData3DDeformation<float> MotionProcessor::get_interpolated_deformation_field(const SignalAxisType signal, const bool cyclic) const
{
	if( this->temp_mvf_filenames_.size() == 0 && this->displacement_fields_.size() == 0)
		throw std::runtime_error("Before calling get_interpolated_deformation_field: Please use prep_displacement_fields() if the fields are not kept in memory, or set_displacement_fields() if they are.");

	if (signal > 1.f || signal< 0.f)
		throw std::runtime_error("Please pass a signal in the range of [0,1].");

	size_t const num_motion_fields = keep_motion_fields_in_memory_? this->displacement_fields_.size() : this->temp_mvf_filenames_.size();
	
	// check in which interval the signal lies
	SignalAxisType signal_on_bin_range;

	if( cyclic )
		signal_on_bin_range = num_motion_fields * signal;
	else
		signal_on_bin_range = (num_motion_fields  - 1)* signal;

	int const bin_floor = int( signal_on_bin_range + 1) -1;
	int const bin_ceil  = int( signal_on_bin_range + 1) % num_motion_fields;

	SignalAxisType const linear_interpolation_weight = signal_on_bin_range - bin_floor;

	/// Constructor
    sirf::ImageWeightedMean<float> dvf_interpolator;

    if(keep_motion_fields_in_memory_)
	{
		dvf_interpolator.add_image( this->displacement_fields_[bin_floor], 1 - linear_interpolation_weight);
	    dvf_interpolator.add_image( this->displacement_fields_[bin_ceil], linear_interpolation_weight);
	}
	else 
	{
		dvf_interpolator.add_image( temp_mvf_filenames_[bin_floor], 1 - linear_interpolation_weight);
	    dvf_interpolator.add_image( temp_mvf_filenames_[bin_ceil], linear_interpolation_weight);
	} 


    dvf_interpolator.process();

    sirf::NiftiImageData3DDisplacement<float> interpolated_dvf( *dvf_interpolator.get_output_sptr() );

    return interpolated_dvf.get_as_deformation_field( interpolated_dvf );
}





std::string MotionProcessor::setup_tmp_folder_name()
{
	std::string const current_folder_prefix = "temp_folder_motion_dyn_";
	std::stringstream tmp_stream;
	tmp_stream << this->temp_folder_prefix_ << current_folder_prefix << this->which_motion_processor_am_i_;
	return tmp_stream.str();

}

std::string MotionProcessor::setup_gt_folder_name()
{
	std::string const gt_folder_prefix = "ground_truth_folder_motion_dyn_";
	std::stringstream name_stream;
	name_stream << this->temp_folder_prefix_ << gt_folder_prefix << this->which_motion_processor_am_i_;
	return name_stream.str();
}


bool MotionProcessor::make_ground_truth_folder() const
{
	try
	{
		std::cout << "Generating ground truth folder " << this->ground_truth_folder_name_ << std::endl;
		boost::filesystem::path dir_to_make(this->ground_truth_folder_name_.c_str());
		bool folder_creation_worked = boost::filesystem::create_directories(dir_to_make);

		return folder_creation_worked;

	}
	catch(boost::system::error_code& e)
	{
		std::cout << e.message() << std::endl;
		throw e;	
	}
}


bool MotionProcessor::make_temp_folder() const
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

bool MotionProcessor::delete_temp_folder() const
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


void MotionProcessor::set_displacement_fields(const std::vector< sirf::NiftiImageData3DDisplacement <float> > &input_displacement_fields)
{
	for(size_t i=0; i<input_displacement_fields.size(); i++)
		this->displacement_fields_.push_back( input_displacement_fields[i] );

	// if( false ) 
	// {
	// 	std::cout << "WARNING: Emtpying the passed vector to save memory." <<std::endl;
	// 	std::vector< sirf::NiftiImageData3DDisplacement<float> >  empty_container;
	// 	input_displacement_fields.swap(empty_container);
	// }
}


sirf::NiftiImageData3DDisplacement<float> MotionProcessor::scale_displacementfields_to_mm( const sirf::NiftiImageData3DDisplacement<float> &dvf ) const
{
    const int* dvf_dims = dvf.get_dimensions() ;

    sirf::VoxelisedGeometricalInfo3D::Spacing voxel_sizes = dvf.get_geom_info_sptr()->get_spacing();

    if( dvf_dims[5] != 3)
     	throw std::runtime_error( "The dimensions of your dvf are not 3D in the 4th spot of the dims but instead." );

    sirf::NiftiImageData3DDisplacement<float> scaled_dvf(dvf);

    for(int nz=0; nz<(int)dvf_dims[3]; nz++)
	for(int ny=0; ny<(int)dvf_dims[2]; ny++)
	for(int nx=0; nx<(int)dvf_dims[1]; nx++)
	{
		
		for(int nv=0; nv<dvf_dims[5]; nv++)
		{
			const int idx[7] = {nx, ny, nz, 0, nv, 0};
			scaled_dvf(idx) =  voxel_sizes[nv] * dvf(idx);
			// scaled_dvf(idx) =  1 * dvf(idx);
		}
	}
	
	return scaled_dvf;
}

void MotionProcessor::prep_displacement_fields()
{

	if(this->displacement_fields_.size() == 0)
		throw std::runtime_error("Please call set_displacement_fields() first.");

	std::cout << "Preparing displacement fields ... " <<std::endl;

	if( !keep_motion_fields_in_memory_ )
	{
		bool const temp_folder_creation_successful = this->make_temp_folder();

		if( temp_folder_creation_successful )
		{
			for(int i=0; i<this->displacement_fields_.size(); i++)
			{
				std::stringstream temp_filename_mvf;
				temp_filename_mvf << this->temp_folder_name_ << this->temp_mvf_prefix_ << i; 
		
				this->displacement_fields_[i].write( temp_filename_mvf.str() );

				temp_filename_mvf << ".nii";
				this->temp_mvf_filenames_.push_back(temp_filename_mvf.str());
			}	
		}

		else
			throw std::runtime_error("The parent directory generation failed. Give a path to which thou hast access rights. Or maybe the directory already exists. This is dangerous. Then you should definitely choose a different temporary folder name.");

		std::vector<MotionFieldType> empty_container;
		this->displacement_fields_.swap(empty_container);
	}

	std::cout << "... finished." <<std::endl;
}

void MotionProcessor::save_ground_truth_displacements( const std::vector< SignalAxisType >& gt_signal_points, const bool cyclic) const
{
	bool const correct_for_offset = false;

	this->make_ground_truth_folder();

	auto minimum_pos = std::min_element(std::begin(gt_signal_points), std::end(gt_signal_points));
	SignalAxisType gt_signal_offset = *minimum_pos;

	std::cout << "Subtracting offset of signal " << gt_signal_offset << std::endl;

	NiftiImageData3DDeformation<float> offset_deformation = this->get_interpolated_deformation_field( gt_signal_offset, cyclic); 

	std::shared_ptr<const NiftiImageData3DDeformation<float> > sptr_inverse_offset_deformation; 

	if( correct_for_offset )
		sptr_inverse_offset_deformation = std::make_shared<const NiftiImageData3DDeformation<float> >(calc_inverse_offset_deformation( offset_deformation ));

	for( size_t i=0; i<gt_signal_points.size(); i++)
	{

		NiftiImageData3DDeformation<float> gt_deformation_with_offset = this->get_interpolated_deformation_field( gt_signal_points[i], cyclic);
		
		std::vector< std::shared_ptr<const Transformation<float> > > vec_gt_def_with_offset;

		vec_gt_def_with_offset.push_back(std::make_shared<const NiftiImageData3DDeformation<float> >( gt_deformation_with_offset ));

		if( correct_for_offset )
			vec_gt_def_with_offset.push_back(sptr_inverse_offset_deformation);
		
		NiftiImageData3DDeformation<float> gt_deformation;
		
		if( correct_for_offset )
			gt_deformation = NiftiImageData3DDeformation<float>::compose_single_deformation(vec_gt_def_with_offset, *sptr_inverse_offset_deformation);
		else
			gt_deformation = gt_deformation_with_offset;
						
		stringstream sstream_output;
		sstream_output << this->ground_truth_folder_name_ << "/gt_deformation_state_" << gt_signal_points[i];
		std::cout << sstream_output.str() << std::endl;

		NiftiImageData3DDisplacement<float> const gt_displacement( gt_deformation );
		gt_displacement.write( sstream_output.str() );		

  	}
}



NiftiImageData3DDeformation<float>  
MotionProcessor::calc_inverse_offset_deformation(NiftiImageData3DDeformation<float> offset_deformation) const
{
	std::shared_ptr<nifti_image> sptr_offset_deformation = offset_deformation.get_raw_nifti_sptr();
	std::shared_ptr<nifti_image> sptr_inverse_offset_deformation = std::make_shared< nifti_image >(*sptr_offset_deformation);

	float const inversion_tolerance = 0.01f;

	std::cout << "Inverting offset vectorfield ... " << std::endl;

	reg_defFieldInvert(sptr_offset_deformation.get(), sptr_inverse_offset_deformation.get(), inversion_tolerance );

	std::cout << " ... finished." << std::endl;

	NiftiImageData3DDeformation<float> inverse_offset( *sptr_inverse_offset_deformation );		

	return inverse_offset;
}


void MRMotionDynamic::bin_mr_acquisitions( MRAcquisitionData& all_acquisitions )
{
	std::cout << "Binning motion dynamics\n";
	
	if(bp_.get_bins().size() == 1)
	{
		std::cout << "We have only one bin, we take all the acquisitions" << std::endl;
		AcquisitionsVector av(all_acquisitions.acquisitions_info());
		ISMRMRD::Acquisition acq;
		for(int i=0; i<all_acquisitions.number(); ++i)
		{
			all_acquisitions.get_acquisition(i,acq);
			acq.resize(1,1,0);
			av.append_acquisition(acq);
		}

		binned_mr_acquisitions_.push_back(av);
		return;
	}

	if(sp_.is_empty())
		throw std::runtime_error( "Please set a signal first. Otherwise you cannot bin your data, you dummy!" );

	all_acquisitions.sort_by_time();

	ISMRMRD::Acquisition acq;
	all_acquisitions.get_acquisition(0, acq);

	TimeAxisType time_offset_seconds = SIRF_SCANNER_MS_PER_TIC/1000.f * acq.getHead().acquisition_time_stamp;

	std::deque< size_t > relevant_acq_numbers;
	std::deque< size_t > acq_not_binned;


	for( size_t i=0; i<all_acquisitions.items(); i++)
		relevant_acq_numbers.push_back( i );

	const std::vector<SignalBin> signal_bins = bp_.get_bins();

	for( int i_bin=0; i_bin<signal_bins.size(); i_bin++)
	{

		auto bin = signal_bins[i_bin];
	
		AcquisitionsVector curr_acq_vector;
		curr_acq_vector.copy_acquisitions_info(all_acquisitions);
		
		acq_not_binned.clear();

		while( relevant_acq_numbers.size() > 0 )	
		{
			auto curr_pos = relevant_acq_numbers[0];
			relevant_acq_numbers.pop_front();	
			ISMRMRD::Acquisition acq;
			all_acquisitions.get_acquisition( curr_pos, acq );
			
			TimeAxisType acq_time_seconds = SIRF_SCANNER_MS_PER_TIC/1000.f * (TimeAxisType)acq.getHead().acquisition_time_stamp - time_offset_seconds;
			
			SignalAxisType signal_of_acq = this->interpolate_signal(acq_time_seconds);
			
			if( is_in_bin(signal_of_acq, bin) )
				curr_acq_vector.append_acquisition(acq);
			else
				acq_not_binned.push_back(curr_pos);					
		}
	
		relevant_acq_numbers.swap(acq_not_binned);
		this->binned_mr_acquisitions_.push_back( curr_acq_vector );
	}
}

std::vector<AcquisitionsVector> MRContrastDynamic::binned_mr_acquisitions_ = std::vector< AcquisitionsVector >(0);

void MRContrastDynamic::bin_mr_acquisitions( MRAcquisitionData& all_acquisitions )
{
	std::cout << "Binning contrast dynamics\n";
	if( true ) //this loop just for RAII reasons to free data
	{
		std::vector<AcquisitionsVector> empty_vec;
		this->binned_mr_acquisitions_.swap( empty_vec );

		cp_.empty_timepoints();
	}

	all_acquisitions.sort_by_time(); 

	size_t const num_acquis = all_acquisitions.number();	
	ISMRMRD::Acquisition acq;
	all_acquisitions.get_acquisition(num_acquis-1, acq);
	float const t_fin = acq.getHead().acquisition_time_stamp;

	all_acquisitions.get_acquisition(0, acq);
	float const t_start = acq.getHead().acquisition_time_stamp;

	TimeAxisType total_time_seconds = SIRF_SCANNER_MS_PER_TIC/1000.f * (t_fin - t_start);
							   
	std::vector< size_t > index_lims;
	const std::vector<SignalBin> signal_bins = bp_.get_bins();

	for( size_t i=0; i<signal_bins.size();i++)
	{
		index_lims.push_back( std::get<2>(signal_bins[i]) *num_acquis );
		cp_.add_timepoint(total_time_seconds * std::get<1>(signal_bins[i]));
	}

	int start_index = 0;
	int stop_index  = 0;

	int const num_bins = bp_.get_num_bins();

	for( size_t i_bin=0; i_bin<num_bins; i_bin++)
	{

		stop_index = ( index_lims[i_bin] < num_acquis ) ? index_lims[i_bin] : num_acquis;

		sirf::AcquisitionsVector av;
		av.copy_acquisitions_info(all_acquisitions);

		for(size_t i=start_index; i<stop_index; i++)
		{
			ISMRMRD::Acquisition acq;
			all_acquisitions.get_acquisition( i, acq );
			acq.resize(1,1,0);
			av.append_acquisition( acq );
		}
		
		this->binned_mr_acquisitions_.push_back( av );
		start_index = stop_index;
	}
}

// ++++++++++++++++++++++++++++++++ PET ++++++++++++++++++++++++++++++++++++++++++++++



TimeBin intersect_time_intervals( const TimeBin& one_interval, const TimeBin& other_interval)
{
	return intersect_intervals<TimeAxisType>(one_interval, other_interval);
}

TimeBinSet intersect_time_bin_sets( const TimeBinSet& one_set, const TimeBinSet& other_set)
{
	TimeBinSet intersected_set;
	for(size_t i=0; i<one_set.size();i++ )
	for(size_t j=0; j<other_set.size();j++ )
	{
		TimeBin temp_intersect = intersect_time_intervals(one_set[i], other_set[j]);
		if( !temp_intersect.is_empty() )			
			intersected_set.push_back( temp_intersect );
	}
	return intersected_set;
}


TimeAxisType get_total_time_in_set(TimeBinSet& set_of_bins )
{
	TimeAxisType t=0;
	for(size_t i_bin=0; i_bin<set_of_bins.size(); i_bin++)	
	{
		t += (set_of_bins[i_bin].max_ - set_of_bins[i_bin].min_);
	}
	return t;
}

void PETDynamic::bin_total_time_interval(TimeBin time_interval_total_dynamic_process_seconds)
{
	if(sp_.is_empty())
		throw std::runtime_error( "Please set a signal first. Otherwise you cannot bin your data, you dummy!" );

	// up-sample the time point to temporal resolution of 1ms
	TimeAxisType delta_time_s = 0.001;

	std::vector<TimeAxisType> upsampled_time_pts;
	std::vector<SignalAxisType> upsampled_signal;

	TimeAxisType current_time = time_interval_total_dynamic_process_seconds.min_;
	std::cout << "Upsampling signal to temporal resolution of 1ms." << std::endl;
	
	while(current_time <= time_interval_total_dynamic_process_seconds.max_)
	{
		upsampled_time_pts.push_back(current_time);
		upsampled_signal.push_back( sp_.linear_interpolate_signal(current_time));

		current_time +=delta_time_s;
	}

	std::cout << "Finished upsampling signal." << std::endl;

	size_t const num_bins = bp_.get_num_bins();
	const std::vector<SignalBin> signal_bins = bp_.get_bins();

	for( size_t i_bin=0; i_bin<num_bins; i_bin++)
	{
		std::cout << "Calculating time bins for " << i_bin << "/ " << num_bins <<std::endl;
		SignalBin bin = signal_bins[i_bin];

		TimeBinSet time_intervals_for_bin;

		auto bin_min = std::get<0>(bin);
		auto bin_max = std::get<2>(bin);

		TimeBin curr_time_interval;
		bool new_bin = true;
		bool no_intersection = true;
		TimeAxisType time_in_bin = 0;

		for(size_t j=0; j<upsampled_signal.size(); ++j)
		{
			SignalAxisType curr_sig = upsampled_signal[j];
			
			if( curr_sig>= bin_min && curr_sig< bin_max)
			{
				if(new_bin)
				{
					curr_time_interval.min_ = upsampled_time_pts[j];
					new_bin = false;
				}

				time_in_bin += delta_time_s;
			}
			else
			{
				no_intersection = false;
				curr_time_interval.max_ = curr_time_interval.min_ + time_in_bin;
				time_intervals_for_bin.push_back( curr_time_interval );
				time_in_bin = 0;
				new_bin = true;
			}
		}

		if( no_intersection )
		{
			curr_time_interval.max_ = curr_time_interval.min_ + time_in_bin;
			time_intervals_for_bin.push_back( curr_time_interval );
		}
		
		this->binned_time_intervals_.push_back( time_intervals_for_bin );
	}
}


TimeAxisType get_time_from_between_two_signal_points(SignalAxisType signal, SignalPoint left_point, SignalPoint right_point)
{
	if(std::abs(right_point.second - left_point.second) < 1e-8)
		return (right_point.first + left_point.first)/TimeAxisType(2);
	else
		return (signal-left_point.second) * (right_point.first - left_point.first) / (right_point.second - left_point.second) + left_point.first;
}

TimeBinSet PETDynamic::get_time_bin_set_for_state( unsigned int const which_state )
{
	if(which_state >= binned_time_intervals_.size())
		throw std::runtime_error( " Please give a number not larger than the number of dynamic states-1");


	return this->binned_time_intervals_[which_state];	
}

TimeAxisType PETDynamic::get_time_spent_in_bin(unsigned int const which_state )
{
	if(which_state >= binned_time_intervals_.size())
		throw std::runtime_error( " Please give a number not larger than the number of dynamic states-1");


	return get_total_time_in_set( this->binned_time_intervals_[which_state] );

}


// void MotionProcessor::align_motion_fields_with_image( const sirf::STIRImageData& img )
// {

// 	size_t const num_disp_fields = this->displacement_fields_.size();

// 	if( num_disp_fields ==0 )
// 		throw std::runtime_error("Please call prep_displacement_fields() first.");

// 	NiftiImageData3D<float> sirf_img( img );
// 	auto sptr_pet_nifti = sirf_img.get_raw_nifti_sptr();

// 	float const img_off_x = sptr_pet_nifti->qoffset_x;
// 	float const img_off_y = sptr_pet_nifti->qoffset_y;
// 	float const img_off_z = sptr_pet_nifti->qoffset_z;
	 
// 	float const img_quart_b = sptr_pet_nifti->quatern_b;
// 	float const img_quart_c = sptr_pet_nifti->quatern_c;
// 	float const img_quart_d = sptr_pet_nifti->quatern_d;
// 	float const img_quart_ac = sptr_pet_nifti->qfac;

// 	for(size_t i=0; i<num_disp_fields; i++)
// 	{

// 		auto sptr_mvf_nifti = this->displacement_fields_[i].get_raw_nifti_sptr();

// 		sptr_mvf_nifti->qoffset_x = img_off_x;
// 		sptr_mvf_nifti->qoffset_y = img_off_y;
// 		sptr_mvf_nifti->qoffset_z = img_off_z;

// 		sptr_mvf_nifti->quatern_b = img_quart_b ;
// 		sptr_mvf_nifti->quatern_c = img_quart_c ;
// 		sptr_mvf_nifti->quatern_d = img_quart_d ;
// 		sptr_mvf_nifti->qfac	  = img_quart_ac;

// 		this->displacement_fields_[i] = NiftiImageData3DDisplacement<float>(*sptr_mvf_nifti);
// 	}
// }