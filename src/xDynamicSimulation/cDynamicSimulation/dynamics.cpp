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

#include "sirf/cReg/NiftyResample.h"
#include "sirf/cReg/NiftiImageData3D.h"
#include "sirf/cReg/NiftiImageData3DDeformation.h"
#include "sirf/cReg/NiftiImageData3DDisplacement.h"

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


 
AcquisitionsVector intersect_mr_acquisition_data( AcquisitionsVector& one_dat, AcquisitionsVector& other_dat)
{

	bool one_dat_is_smaller = ( one_dat.items() >= other_dat.items() );

	typedef std::vector<uint32_t> CounterBox;

	CounterBox one_counters, other_counters;


	for( size_t i=0; i<one_dat.items(); i++)
	{
		auto sptr_acq = one_dat.get_acquisition_sptr( i );
		one_counters.push_back(sptr_acq->getHead().scan_counter);
	}

	for( size_t i=0; i<other_dat.items(); i++)
	{
		ISMRMRD::Acquisition acq;

		auto sptr_acq = other_dat.get_acquisition_sptr( i );
		other_counters.push_back(sptr_acq->getHead().scan_counter);
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
		auto sptr_acq = smaller_data_container.get_acquisition_sptr( i );
		uint32_t acquis_counter = sptr_acq->getHead().scan_counter;
		if(std::find(intersected_counters.begin(), intersected_counters.end(), acquis_counter) != intersected_counters.end()) 
		{
			intersection.append_acquisition_sptr(sptr_acq);
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
	std::cout << "size in the getter " << epiph( this->binned_mr_acquisitions_.size()) <<std::endl;
	return this->binned_mr_acquisitions_;
};

sirf::AcquisitionsVector aMRDynamic::get_binned_mr_acquisitions( int const bin_num )
{
	if(bin_num >= this->num_simul_states_)
		throw std::runtime_error("Please access only bin numbers in the range of 0 and num_simul_states_-1.");
	
	return this->binned_mr_acquisitions_[bin_num];
};


// static member variable initialization
int ContrastDynamic::num_simul_states_ = 0;
std::vector< TimeAxisType > ContrastDynamic::time_points_sampled_ = std::vector<TimeAxisType>(0);


ContrastDynamic::ContrastDynamic(int const num_simul_states): aDynamic()
{
	this->num_simul_states_ = num_simul_states;
	this->set_bins( num_simul_states );
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

void ContrastDynamic::set_bins( int const num_bins )
{
	if(num_simul_states_ != 1)
		this->num_simul_states_ = num_bins+1;

	for(int i_state=0; i_state<this->num_simul_states_; i_state++)
	{	
		SignalBin bin;

		std::get<0>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins)- 1.f/(2*num_bins);;
		std::get<1>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins);
		std::get<2>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins)+ 1.f/(2*num_bins);
	
		this->signal_bins_.push_back( bin );
	}
}


int MotionDynamic::num_total_motion_dynamics_ = 0;

MotionDynamic::MotionDynamic():aDynamic()
{
	this->which_motion_dynamic_am_i_ = num_total_motion_dynamics_;
	this->num_total_motion_dynamics_ += 1;
	
	this->temp_folder_name_ = setup_tmp_folder_name();
	this->ground_truth_folder_name_ = setup_gt_folder_name();
}

MotionDynamic::MotionDynamic(int const num_simul_states) : aDynamic()
{
	this->num_simul_states_ =num_simul_states;
	this->set_bins(num_simul_states_);

	this->which_motion_dynamic_am_i_ = num_total_motion_dynamics_;
	this->num_total_motion_dynamics_ += 1;

	this->temp_folder_name_ = setup_tmp_folder_name();
	this->ground_truth_folder_name_ = setup_gt_folder_name();
}


MotionDynamic::~MotionDynamic()
{ 
	// if( this->destroy_upon_deletion_)
	// this->delete_temp_folder();

	this->num_total_motion_dynamics_ -= 1; 
}


NiftiImageData3DDeformation<float> MotionDynamic::get_interpolated_deformation_field(SignalAxisType signal)
{
	if( this->temp_mvf_filenames_.size() == 0 && this->displacement_fields_.size() == 0)
		throw std::runtime_error("Before calling get_interpolated_deformation_field: Please use prep_displacement_fields() if the fields are not kept in memory, or set_displacement_fields() if they are.");

	if (signal > 1.f || signal< 0.f)
		throw std::runtime_error("Please pass a signal in the range of [0,1].");

	size_t const num_motion_fields = keep_motion_fields_in_memory_? this->displacement_fields_.size() : this->temp_mvf_filenames_.size();
	
	// check in which interval the signal lies
	SignalAxisType signal_on_bin_range;

	if( this->is_cyclic_dynamic_ )
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




int MotionDynamic::get_which_motion_dynamic_am_i(){ return this->which_motion_dynamic_am_i_; }
int MotionDynamic::get_num_total_motion_dynamics(){ return this->num_total_motion_dynamics_; }

std::string MotionDynamic::setup_tmp_folder_name()
{
	std::string const current_folder_prefix = "temp_folder_motion_dyn_";
	std::stringstream tmp_stream;
	tmp_stream << this->temp_folder_prefix_ << current_folder_prefix << this->which_motion_dynamic_am_i_;
	return tmp_stream.str();

}

std::string MotionDynamic::setup_gt_folder_name()
{
	std::string const gt_folder_prefix = "ground_truth_folder_motion_dyn_";
	std::stringstream name_stream;
	name_stream << this->temp_folder_prefix_ << gt_folder_prefix << this->which_motion_dynamic_am_i_;
	return name_stream.str();
}

std::string MotionDynamic::get_temp_folder_name()
{
	return this->temp_folder_name_;
}

bool MotionDynamic::make_ground_truth_folder()
{
	try
	{
		std::cout << "Generating temporary folder " << this->ground_truth_folder_name_ << std::endl;
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

// delete this!
// void MotionDynamic::set_displacement_fields( ISMRMRD::NDArray< DataTypeMotionFields >& motion_fields, bool const motion_fields_are_cyclic)
// {
	
// 	if ( motion_fields_are_cyclic )
// 	{
// 		this->is_cyclic_dynamic_ = true;
// 		this->set_bins( this->num_simul_states_ );
// 	}	


// 	using namespace ISMRMRD;

// 	const size_t* dimensions = motion_fields.getDims();

// 	size_t const Nt = dimensions[0];
// 	size_t const Nv = dimensions[1];
// 	size_t const Nz = dimensions[2];
// 	size_t const Ny = dimensions[3];
// 	size_t const Nx = dimensions[4];

// 	for(size_t nt=0; nt<Nt; nt++)
// 	{
		
// 		Image<DataTypeMotionFields> img(dimensions[4],dimensions[3], dimensions[2], dimensions[1]);
 		
//  		for(uint16_t  nv= 0; nv<Nv ; nv++)
// 		for(uint16_t  nz= 0; nz<Nz ; nz++)
// 		for(uint16_t  ny= 0; ny<Ny ; ny++)
// 		for(uint16_t  nx= 0; nx<Nx ; nx++)
// 		{
// 			// size_t const lin_index = ((((Nt-1 -nt)*Nv + Nv-1 -nv)*Nz + Nz-1 - nz)*Ny + Ny-1 - ny)*Nx + Nx-1 - nx;
// 			size_t const lin_index = (((nt*Nv + nv)*Nz + nz)*Ny + ny)*Nx + nx;
// 			img(nx,ny,nz,nv) = 	  *(motion_fields.begin() + lin_index);
// 		}
// 		// this->displacement_fields_.push_back(img);
// 	}
// }


void MotionDynamic::set_ground_truth_folder_name( std::string const name_existing_folder_path )
{
	this->ground_truth_folder_name_ = name_existing_folder_path;
	this->make_ground_truth_folder();
}

void MotionDynamic::set_displacement_fields( std::vector< sirf::NiftiImageData3DDisplacement <float> > &input_displacement_fields, bool const motion_fields_are_cyclic)
{
	if ( motion_fields_are_cyclic )
	{
		this->is_cyclic_dynamic_ = true;
		this->set_bins( this->num_simul_states_ );
	}

	for(size_t i=0; i<input_displacement_fields.size(); i++)
		this->displacement_fields_.push_back( input_displacement_fields[i] );

	if( true ) 
	{
		std::cout << "WARNING: Emtpying the passed vector to save memory." <<std::endl;
		std::vector< sirf::NiftiImageData3DDisplacement<float> >  empty_container;
		input_displacement_fields.swap(empty_container);
	}

}

sirf::NiftiImageData3DDisplacement<float> MotionDynamic::scale_displacementfields_to_mm( const sirf::NiftiImageData3DDisplacement<float> &dvf )
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

void MotionDynamic::prep_displacement_fields()
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

void MotionDynamic::save_ground_truth_displacements( std::vector< SignalAxisType > gt_signal_points)
{
	bool const correct_for_offset = false;


	this->make_ground_truth_folder();

	std::vector<SignalAxisType>::iterator minimum_pos = std::min_element(std::begin(gt_signal_points), std::end(gt_signal_points));
	SignalAxisType gt_signal_offset = *minimum_pos;

	std::cout << "Subtracting offset of signal " << gt_signal_offset << std::endl;

	NiftiImageData3DDeformation<float> offset_deformation = this->get_interpolated_deformation_field( gt_signal_offset ); 


	std::shared_ptr<const NiftiImageData3DDeformation<float> > sptr_inverse_offset_deformation; 
	//= std::make_shared<const NiftiImageData3DDeformation<float> >(calc_inverse_offset_deformation( offset_deformation ));
	if( correct_for_offset )
		sptr_inverse_offset_deformation = std::make_shared<const NiftiImageData3DDeformation<float> >(calc_inverse_offset_deformation( offset_deformation ));

	for( size_t i=0; i<gt_signal_points.size(); i++)
	{

		NiftiImageData3DDeformation<float> gt_deformation_with_offset = this->get_interpolated_deformation_field( gt_signal_points[i] );
		
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

		// NiftiImageData3DDisplacement<float> const gt_displacements_with_offset( gt_deformation_with_offset );
		// gt_displacements_with_offset.write( sstream_output.str() );		
  	}
}

void MotionDynamic::save_ground_truth_displacements( void )
{
	std::vector< SignalBin > simulated_motion_bins = this->get_bins();
	std::vector< SignalAxisType > bin_centers;

	for( size_t i=0; i<simulated_motion_bins.size(); i++)
	{	
		std::cout << "pushing center " << std::get<1>( simulated_motion_bins[i]) << std::endl;
		bin_centers.push_back( std::get<1>( simulated_motion_bins[i]) );
	}

	this->save_ground_truth_displacements( bin_centers );

}

NiftiImageData3DDeformation<float>  
MotionDynamic::calc_inverse_offset_deformation( NiftiImageData3DDeformation<float> offset_deformation )
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


void MRMotionDynamic::bin_mr_acquisitions( AcquisitionsVector& all_acquisitions )
{

	if(this->dyn_signal_.size() == 0)
		throw std::runtime_error( "Please set a signal first. Otherwise you cannot bin your data, you dummy!" );

	AcquisitionsVector time_ordered_acquisitions = all_acquisitions;
	time_ordered_acquisitions.time_order();
	TimeAxisType time_offset = SIRF_SCANNER_MS_PER_TIC * time_ordered_acquisitions.get_acquisition_sptr(0)->acquisition_time_stamp();


	std::deque< size_t > relevant_acq_numbers;
	std::deque< size_t > acq_not_binned;


	for( size_t i=0; i<all_acquisitions.items(); i++)
		relevant_acq_numbers.push_back( i );

	for( int i_bin=0; i_bin<this->signal_bins_.size(); i_bin++)
	{

		auto bin = this->signal_bins_[i_bin];
	
		AcquisitionsVector curr_acq_vector;
		curr_acq_vector.copy_acquisitions_info(all_acquisitions);
		
		acq_not_binned.clear();

		while( relevant_acq_numbers.size() > 0 )	
		{
			auto curr_pos = relevant_acq_numbers[0];
			relevant_acq_numbers.pop_front();	
			
			auto sptr_acq = all_acquisitions.get_acquisition_sptr( curr_pos );
			auto acq_hdr = sptr_acq->getHead();
			
			TimeAxisType acq_time = SIRF_SCANNER_MS_PER_TIC * (TimeAxisType)acq_hdr.acquisition_time_stamp - time_offset;
			
			SignalAxisType signal_of_acq = this->linear_interpolate_signal( acq_time );
			
			if( is_in_bin(signal_of_acq, bin) )
				curr_acq_vector.append_acquisition_sptr(sptr_acq);
			else
				acq_not_binned.push_back(curr_pos);					
		}
	
		relevant_acq_numbers.swap(acq_not_binned);
		this->binned_mr_acquisitions_.push_back( curr_acq_vector );
	}
}


// void MRMotionDynamic::prep_displacement_fields()
// {
// 	MotionDynamic::prep_displacement_fields();
// }

std::vector<AcquisitionsVector> MRContrastDynamic::binned_mr_acquisitions_ = std::vector< AcquisitionsVector >(0);

std::vector<sirf::AcquisitionsVector> MRContrastDynamic::get_binned_mr_acquisitions( void )
{
	std::cout << "size in the getter " << epiph( this->binned_mr_acquisitions_.size()) <<std::endl;
	return this->binned_mr_acquisitions_;
};

sirf::AcquisitionsVector MRContrastDynamic::get_binned_mr_acquisitions( int const bin_num )
{
	if(bin_num >= this->num_simul_states_)
		throw std::runtime_error("Please access only bin numbers in the range of 0 and num_simul_states_-1.");
	
	return this->binned_mr_acquisitions_[bin_num];
};

void MRContrastDynamic::bin_mr_acquisitions( AcquisitionsVector& all_acquisitions )
{
	if( true ) // empty the old bins -> every mr contrast dynamic must hold same binned data
	{
		std::vector<AcquisitionsVector> empty_vec;
		this->binned_mr_acquisitions_.swap( empty_vec );

		this->time_points_sampled_.clear();
	}

	AcquisitionsVector time_ordered_acquisitions = all_acquisitions;
	time_ordered_acquisitions.time_order(); 

	size_t const num_acquis = time_ordered_acquisitions.number();	

	TimeAxisType total_time = SIRF_SCANNER_MS_PER_TIC * (time_ordered_acquisitions.get_acquisition_sptr(num_acquis-1)->acquisition_time_stamp() -
							   time_ordered_acquisitions.get_acquisition_sptr(0)->acquisition_time_stamp());

	std::vector< size_t > index_lims;

	for( size_t i=0; i<this->num_simul_states_;i++)
	{
		index_lims.push_back( std::get<2>(this->signal_bins_[i]) *num_acquis );
		this->time_points_sampled_.push_back( total_time * std::get<1>(this->signal_bins_[i]) );
	}

	int start_index = 0;
	int stop_index  = 0;

	for( size_t i_bin=0; i_bin<this->num_simul_states_; i_bin++)
	{

		stop_index = ( index_lims[i_bin] < num_acquis ) ? index_lims[i_bin] : num_acquis;

		sirf::AcquisitionsVector av;
		av.copy_acquisitions_info(time_ordered_acquisitions);

		for(size_t i=start_index; i<stop_index; i++)
		{
			auto sptr_acq = time_ordered_acquisitions.get_acquisition_sptr( i );
			av.append_acquisition_sptr( sptr_acq );
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

aPETDynamic::aPETDynamic(int const num_simul_states): aDynamic(num_simul_states){}


void aPETDynamic::bin_total_time_interval(TimeBin time_interval_total_dynamic_process)
{
	if(this->dyn_signal_.size() == 0)
		throw std::runtime_error( "Please set a signal first. Otherwise you cannot bin your data, you dummy!" );

	// up-sample the time point to temporal resolution of 1ms
	TimeAxisType delta_time_ms = 1;

	std::vector<TimeAxisType> upsampled_time_pts;
	std::vector<SignalAxisType> upsampled_signal;

	TimeAxisType current_time = time_interval_total_dynamic_process.min_;
	std::cout << "Upsampling signal to temporal resolution of 1ms." << std::endl;
	
	while(current_time <= time_interval_total_dynamic_process.max_)
	{
		upsampled_time_pts.push_back(current_time);
		upsampled_signal.push_back( this->linear_interpolate_signal(current_time));

		current_time +=delta_time_ms;
	}

	std::cout << "Finished upsampling signal." << std::endl;

	size_t const num_bins = signal_bins_.size();

	for( size_t i_bin=0; i_bin<num_bins; i_bin++)
	{
		std::cout << "Calculating time bins for " << i_bin << "/ " << num_bins <<std::endl;
		SignalBin bin = this->signal_bins_[i_bin];

		TimeBinSet time_intervals_for_bin;

		auto bin_min = std::get<0>(bin);
		auto bin_max = std::get<2>(bin);

		TimeBin curr_time_interval;
		bool new_bin = true;
		TimeAxisType time_in_bin = 0;
		for(size_t j=0; j<upsampled_signal.size(); ++j)
		{
			SignalAxisType curr_sig = upsampled_signal[j];
			
			if( curr_sig>= bin_min && curr_sig< bin_max)
			{
				if(new_bin)
					curr_time_interval.min_ = upsampled_time_pts[j];

				time_in_bin += delta_time_ms;
			}
			else
			{
				curr_time_interval.max_ = curr_time_interval.min_ + time_in_bin;
				time_intervals_for_bin.push_back( curr_time_interval );
				time_in_bin = 0;
				new_bin = true;
			}
		}
		this->binned_time_intervals_.push_back( time_intervals_for_bin );
		
	}

	/*	
	size_t const num_bins = signal_bins_.size();
	size_t const num_signal_supports = dyn_signal_.size();
	
	TimeAxisType leftmost_left_edge = std::min<TimeAxisType>( time_interval_total_dynamic_process.min_, dyn_signal_[0].first );
	TimeAxisType rightmost_right_edge = std::max<TimeAxisType>( time_interval_total_dynamic_process.max_,  dyn_signal_[num_signal_supports-1].first );

	TimeBinSet temp_set;
	temp_set.push_back(time_interval_total_dynamic_process);
	for( size_t i_bin=0; i_bin<num_bins; i_bin++)
	{
		SignalBin bin = this->signal_bins_[i_bin];

		auto bin_min = std::get<0>(bin);
		auto bin_max = std::get<2>(bin);

		TimeBinSet time_intervals_for_bin;

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

		size_t const num_left_edges = left_bin_edges.size();
		size_t const num_right_edges = right_bin_edges.size();

		if (  std::abs( (float)num_left_edges - (float)num_right_edges )  > 1 )
			throw std::runtime_error( "You got yourself a very weird combination of signal and bin number in your simulation. Consider passing another dynamic signal please.");	
		
		if( num_left_edges ==0 && num_right_edges == 0)
		{
			TimeBin curr_bin(leftmost_left_edge, rightmost_right_edge);
			time_intervals_for_bin.push_back(curr_bin);
		}
		else if( num_left_edges == num_right_edges )
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

			TimeBin leftmost_bin( leftmost_left_edge, right_bin_edges[0]);
			time_intervals_for_bin.push_back(leftmost_bin);

		}
		time_intervals_for_bin = intersect_time_bin_sets(time_intervals_for_bin, temp_set);	
		this->binned_time_intervals_.push_back(time_intervals_for_bin);
	}
	*/

}


TimeAxisType get_time_from_between_two_signal_points(SignalAxisType signal, SignalPoint left_point, SignalPoint right_point)
{
	if(std::abs(right_point.second - left_point.second) < 1e-8)
		return (right_point.first + left_point.first)/TimeAxisType(2);
	else
		return (signal-left_point.second) * (right_point.first - left_point.first) / (right_point.second - left_point.second) + left_point.first;
}

TimeBinSet aPETDynamic::get_time_bin_set_for_state( unsigned int const which_state )
{
	if(which_state >= binned_time_intervals_.size())
		throw std::runtime_error( " Please give a number not larger than the number of dynamic states-1");


	return this->binned_time_intervals_[which_state];	
}

TimeAxisType aPETDynamic::get_time_spent_in_bin(unsigned int const which_state )
{
	if(which_state >= binned_time_intervals_.size())
		throw std::runtime_error( " Please give a number not larger than the number of dynamic states-1");


	return get_total_time_in_set( this->binned_time_intervals_[which_state] );

}




// void PETMotionDynamic::prep_displacement_fields( void )
// {
// 	if(this->displacement_fields_.size() == 0)
// 		throw std::runtime_error("Please call set_displacements_fields() first.");

// 	std::cout << "Preparing PET displacement fields..." <<std::endl;

// 	bool const temp_folder_creation_successful = this->make_temp_folder();

// 	if( temp_folder_creation_successful )
// 	{
// 		for(int i=0; i<this->displacement_fields_.size(); i++)
// 		{
// 			std::stringstream temp_filename_mvf;
// 			temp_filename_mvf << this->get_temp_folder_name() << this->temp_mvf_prefix_ << i;

// 			data_io::write_MVF_from_ISMRMRD_Image_to_Analyze_In_PET_Geometry<DataTypeMotionFields> (temp_filename_mvf.str(), this->displacement_fields_[i]);
// 			temp_filename_mvf << ".hdr";
// 			this->temp_mvf_filenames_.push_back(temp_filename_mvf.str());
// 		}

// 		std::vector<MotionFieldType> empty_container;
// 		this->displacement_fields_.swap(empty_container); 

// 	}
// 	else
// 		throw std::runtime_error("The parent directory generation failed. Give a path to which thou hast access rights. Or maybe the directory already exists. This is dangerous. Then you should definitely choose a different temporary folder name.");

// 	for( size_t i=0; i<temp_mvf_filenames_.size(); i++)
// 	{
// 		NiftiImageData3DDeformation<float> temp_deformation( this->temp_mvf_filenames_[i] );
// 		this->sirf_displacement_fields_.push_back( temp_deformation );
// 	}

// 	this->delete_temp_folder();
// 	std::cout << "... finished." <<std::endl;
// }

void PETMotionDynamic::align_motion_fields_with_image( const sirf::STIRImageData& img )
{

	size_t const num_disp_fields = this->displacement_fields_.size();

	if( num_disp_fields ==0 )
		throw std::runtime_error("Please call prep_displacement_fields() first.");

	NiftiImageData3D<float> sirf_img( img );
	auto sptr_pet_nifti = sirf_img.get_raw_nifti_sptr();

	float const img_off_x = sptr_pet_nifti->qoffset_x;
	float const img_off_y = sptr_pet_nifti->qoffset_y;
	float const img_off_z = sptr_pet_nifti->qoffset_z;
	 
	float const img_quart_b = sptr_pet_nifti->quatern_b;
	float const img_quart_c = sptr_pet_nifti->quatern_c;
	float const img_quart_d = sptr_pet_nifti->quatern_d;
	float const img_quart_ac = sptr_pet_nifti->qfac;


	// float const img_dx = sptr_pet_nifti->dx;
 //    float const img_dy = sptr_pet_nifti->dy;
 //    float const img_dz = sptr_pet_nifti->dz;
 //    float const img_dt = sptr_pet_nifti->dt;
 //    float const img_du = sptr_pet_nifti->du;
 //    float const img_dv = sptr_pet_nifti->dv;
 //    float const img_dw = sptr_pet_nifti->dw;

	for(size_t i=0; i<num_disp_fields; i++)
	{

		auto sptr_mvf_nifti = this->displacement_fields_[i].get_raw_nifti_sptr();

		sptr_mvf_nifti->qoffset_x = img_off_x;
		sptr_mvf_nifti->qoffset_y = img_off_y;
		sptr_mvf_nifti->qoffset_z = img_off_z;

		sptr_mvf_nifti->quatern_b = img_quart_b ;
		sptr_mvf_nifti->quatern_c = img_quart_c ;
		sptr_mvf_nifti->quatern_d = img_quart_d ;
		sptr_mvf_nifti->qfac	  = img_quart_ac;

		this->displacement_fields_[i] = NiftiImageData3DDisplacement<float>(*sptr_mvf_nifti);

	}
}





















