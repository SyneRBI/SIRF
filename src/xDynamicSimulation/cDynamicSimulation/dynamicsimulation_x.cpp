/*
author: johannes mayer
date: 15. March 2018

*/

#include <memory>
#include <mutex>
#include <sstream>

#include "sirf/cDynamicSimulation/dynamicsimulation_x.h"

#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

#include "sirf/cReg/NiftiImageData3DDeformation.h"

#include "sirf/cDynamicSimulation/dynsim_deformer.h"



using namespace sirf;
using std::cout;
using std::endl;


void MRDynamicSimulation::set_filename_rawdata( std::string const filename_template_rawdata ) 
{
	aDynamicSimulation::set_filename_rawdata( filename_template_rawdata );
	this->extract_hdr_information();
}



void MRDynamicSimulation::write_simulation_results( const std::string& filename_output_with_h5_extension ) 
{	
	try	
	{
		cout << "Started writing simulation output to: " << filename_output_with_h5_extension <<endl;
		cout << "Number of acquisitions to write: " << this->target_acquisitions_.number() << endl;

		std::stringstream serialized_hdr;
		ISMRMRD::serialize(this->hdr_, serialized_hdr);
		target_acquisitions_.set_acquisitions_info( serialized_hdr.str() ); 

		target_acquisitions_.time_order();		

		target_acquisitions_.write( filename_output_with_h5_extension.c_str() );
		cout << "Finished writing simulation output."<<endl;
	}
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		cout << "Maybe you forgot to give the correct filename ending .h5 in this case." << endl;
	}
}

void MRDynamicSimulation::save_ground_truth_displacements( void )
{
	for(size_t i=0; i<this->motion_dynamics_.size(); i++)
	{
		this->motion_dynamics_[i]->save_ground_truth_displacements();
	}
}

void MRDynamicSimulation::add_dynamic( std::shared_ptr<MRMotionDynamic> sptr_motion_dyn )
{
	this->motion_dynamics_.push_back( sptr_motion_dyn );
};


void MRDynamicSimulation::add_dynamic( std::shared_ptr<MRContrastDynamic> sptr_contrast_dyn) 
{
	this->contrast_dynamics_.push_back( sptr_contrast_dyn );
};


void MRDynamicSimulation::simulate_dynamics( void )
{

	cout << "Simulating dynamic data acquisition... " <<endl;
	this->simulate_simultaneous_motion_contrast_dynamics();		
	
}

void MRDynamicSimulation::simulate_simultaneous_motion_contrast_dynamics()
{
	cout << "Simulating motion and contrast dynamics... " <<endl;

	this->extract_hdr_information();
	this->acq_model_.setTraj( this->sptr_trajectory_ );

	this->mr_cont_gen_.map_contrast();
		
	// all contrast dynamic variations are sampled at the same timepoint.
    size_t const num_contrast_dyns = this->contrast_dynamics_.size();
    size_t num_contrast_states = (num_contrast_dyns > 0)? contrast_dynamics_[0]->get_num_simul_states() : 1;

    std::vector< int > num_states_per_motion;
    size_t const num_motion_dyns = this->motion_dynamics_.size();
	for(size_t i=0; i<num_motion_dyns; i++)	
	{
		num_states_per_motion.push_back(motion_dynamics_[i]->get_num_simul_states());			
		motion_dynamics_[i]->prep_displacement_fields();
	}

	LinearCombiGenerator lcg(num_states_per_motion);
	
	size_t const num_tot_motion_states = lcg.get_num_total_combinations();

	std::vector< DimensionsType >  all_dynamic_state_combos = lcg.get_all_combinations();

    auto binned_contrast_acquisitions = (num_contrast_dyns > 0 ) ? contrast_dynamics_[0]->get_binned_mr_acquisitions()  : std::vector< AcquisitionsVector >(0);
    auto sampled_contrast_timepoints = (num_contrast_dyns > 0 ) ? contrast_dynamics_[0]->get_time_points_sampled() : std::vector< TimeAxisType >(0);
    
    for( size_t i_dyn_state=0; i_dyn_state < num_tot_motion_states; i_dyn_state++)
	{
		cout << " ++++++++++++    Acquisition dynamic motion state #" << i_dyn_state + 1 << "/" << num_tot_motion_states << "++++++++++++++" << endl;

		DimensionsType current_combination = all_dynamic_state_combos[i_dyn_state];
		AcquisitionsVector acquisitions_for_this_motion_state = this->all_source_acquisitions_;

        for( int i_motion_dyn = 0; i_motion_dyn<num_motion_dyns; i_motion_dyn++ )
		{
			auto motion_dyn = this->motion_dynamics_[ i_motion_dyn ];
			int const current_bin = current_combination[ i_motion_dyn ];						

			AcquisitionsVector acquis_in_motion_bin = motion_dyn->get_binned_mr_acquisitions( current_bin );
			acquisitions_for_this_motion_state = intersect_mr_acquisition_data(acquisitions_for_this_motion_state, acquis_in_motion_bin);

		}
    
        cout << "# of mr acquis in this motion state: " << acquisitions_for_this_motion_state.number() << endl;
    
        if( acquisitions_for_this_motion_state.number() > 0)
        {
            
            for( size_t i_contrast_state=0; i_contrast_state<num_contrast_states; i_contrast_state++)
            {
           		cout << " ++++++++++++    Acquisition dynamic contrast state #" << i_contrast_state + 1 << "/" << num_contrast_states << "++++++++++++++" << endl;
				auto acquisitions_for_this_contrast_state = acquisitions_for_this_motion_state;

				if( binned_contrast_acquisitions.size() >0)
				{
	                AcquisitionsVector acquis_in_contrast_state = binned_contrast_acquisitions[ i_contrast_state ];
	                acquisitions_for_this_contrast_state = intersect_mr_acquisition_data(acquisitions_for_this_contrast_state, acquis_in_contrast_state);           
                }              

                cout << "# of mr acquis in this dynamic motion x contrast state: " << acquisitions_for_this_contrast_state.number() << endl;

                std::vector<sirf::NiftiImageData3DDeformation<float> > all_motion_fields;
                
                for( int i_motion_dyn = 0; i_motion_dyn<num_motion_dyns; i_motion_dyn++ )
                {
                    cout << "Preparing motion fields for motion dynamic # " << i_motion_dyn << endl;

                    auto motion_dyn = this->motion_dynamics_[i_motion_dyn];
                    std::vector< SignalBin > signal_bins = motion_dyn->get_bins();

                    int const motion_bin_number = current_combination[ i_motion_dyn ];						
                    SignalBin bin = signal_bins[ motion_bin_number ];
	                
                    all_motion_fields.push_back( motion_dyn->get_interpolated_deformation_field( std::get<1>(bin) ) ); 

                }

                if( acquisitions_for_this_contrast_state.number() > 0)
                {
                    TimeAxisType current_time_point =  sampled_contrast_timepoints.size()>0 ? sampled_contrast_timepoints[i_contrast_state] : (TimeAxisType)0;
                    
                    for( int i_contrast_dyn = 0; i_contrast_dyn<num_contrast_dyns; i_contrast_dyn++ )
                    {
                        auto sptr_contrast_dyn = this->contrast_dynamics_[i_contrast_dyn];
                        auto contrast_signal = sptr_contrast_dyn->linear_interpolate_signal(current_time_point);
                        TissueParameterList tissueparameter_list_to_replace = sptr_contrast_dyn->get_interpolated_tissue_params( contrast_signal );

                        for( size_t i_tiss=0; i_tiss< tissueparameter_list_to_replace.size(); i_tiss++ )
                        {
                            TissueParameter curr_param = tissueparameter_list_to_replace[i_tiss];
                            this->mr_cont_gen_.replace_petmr_tissue_parameters( curr_param.label_, curr_param );	
                        }
                    }

                    this->mr_cont_gen_.map_contrast();//crucial call, as the deformation results in deformed contrast generator data
                    
                    if( all_motion_fields.size() > 0 )
                        DynamicSimulationDeformer::deform_contrast_generator(this->mr_cont_gen_, all_motion_fields);
                    
                    this->source_acquisitions_ = acquisitions_for_this_contrast_state;
                    this->acquire_raw_data();	
                }
            }
        }
	}
	
    this->noise_generator_.add_noise(this->target_acquisitions_);

	for(size_t i=0; i<num_motion_dyns; i++)
		this->motion_dynamics_[i]->delete_temp_folder();	

}


void MRDynamicSimulation::extract_hdr_information( void )
{
	this->hdr_ = mr_io::read_ismrmrd_header( filename_rawdata_ );

	this->sptr_trajectory_->set_header( this->hdr_ );
	this->sptr_trajectory_->compute_trajectory();

	sptr_trajectory_->overwrite_ismrmrd_trajectory_info( this->hdr_ );

	this->acq_model_.setISMRMRDHeader( this->hdr_ );
	this->mr_cont_gen_.set_rawdata_header( this->hdr_ );
	
}


void MRDynamicSimulation::set_trajectory( std::shared_ptr<aTrajectoryContainer> sptr_trajectory)
{
	this->sptr_trajectory_ = sptr_trajectory;
	sptr_trajectory_->overwrite_ismrmrd_trajectory_info( this->hdr_ );
	
	this->mr_cont_gen_.set_rawdata_header( this->hdr_ );

	this->set_noise_scaling( this->sptr_trajectory_ );

}


void MRDynamicSimulation::set_noise_scaling( std::shared_ptr<aTrajectoryContainer> sptr_trajectory )
{
	ISMRMRD::TrajectoryType traj_type = sptr_trajectory->get_traj_type();

	if( traj_type == ISMRMRD::TrajectoryType::OTHER )
		this->noise_generator_.set_sampling_specific_scaling(RPE_NOISE_SCALING);
}


void MRDynamicSimulation::set_coilmaps( ISMRMRD::Image< complex_float_t > &coilmaps )
{
	this->coilmaps_ = coilmaps;
}

void MRDynamicSimulation::shift_time_start_to_zero( void )
{
	this->all_source_acquisitions_.time_order();
	auto sptr_acquis = all_source_acquisitions_.get_acquisition_sptr(0);
	uint32_t const t0 = sptr_acquis->acquisition_time_stamp();
	for(size_t i=0; i<all_source_acquisitions_.number(); i++)
	{
		sptr_acquis = all_source_acquisitions_.get_acquisition_sptr(i);
		sptr_acquis->acquisition_time_stamp() -= t0;
	}
}

void MRDynamicSimulation::set_all_source_acquisitions(MRDataType& acquisitions )
{
	this->all_source_acquisitions_ = acquisitions;
	this->target_acquisitions_.copy_acquisitions_info( this->all_source_acquisitions_ );

	this->shift_time_start_to_zero();
}


void MRDynamicSimulation::set_SNR(float const SNR)
{
	this->noise_generator_.set_SNR(SNR);
}

void MRDynamicSimulation::set_noise_label(size_t const label)
{
	auto const signal_in_label = this->mr_cont_gen_.get_signal_for_tissuelabel(label);
	auto const abs_signal = std::abs( signal_in_label );
	// float const abs_signal = 1.f;
	std::cout << "Adding signal " << abs_signal << " for label " << label << std::endl;
	this->noise_generator_.set_signal_img( abs_signal );
}


void MRDynamicSimulation::acquire_raw_data( void )
{

	if( this->coilmaps_.getNumberOfDataElements() == 0)
		throw std::runtime_error("Please make sure to set the coilmaps prior to starting the simulation.");

	std::vector< ISMRMRD::Image< complex_float_t> > contrast_filled_volumes = this->mr_cont_gen_.get_contrast_filled_volumes();

	size_t const num_contrasts = contrast_filled_volumes.size();


	size_t Nx = contrast_filled_volumes[0].getMatrixSizeX();
	size_t Ny = contrast_filled_volumes[0].getMatrixSizeY();
	size_t Nz = contrast_filled_volumes[0].getMatrixSizeZ();
	size_t Nc = this->coilmaps_.getNumberOfChannels();

	auto csm = vol_orientator_.reorient_image(this->coilmaps_);
	CoilDataAsCFImage csm_as_img( csm.getMatrixSizeX(), csm.getMatrixSizeY(), csm.getMatrixSizeZ() , Nc);
	csm_as_img.image() = csm;

	unsigned int offset = 0;


	for( size_t i_contrast=0; i_contrast<num_contrasts; i_contrast++)
	{
		cout << "Acquisition of contrast " << i_contrast << endl;
		ISMRMRD::Image<complex_float_t> curr_cont = contrast_filled_volumes[i_contrast];
		curr_cont = vol_orientator_.reorient_image(curr_cont);

		ImageWrap curr_img_wrap(IMG_DATA_TYPE, new ISMRMRD::Image< complex_float_t >(curr_cont));		

		AcquisitionsVector acq_vec;
		acq_vec.copy_acquisitions_info( this->source_acquisitions_ );
		
		for( size_t i_acqu=0; i_acqu<this->source_acquisitions_.items(); i_acqu++)
		{
			auto sptr_acq = this->source_acquisitions_.get_acquisition_sptr(i_acqu);
			ISMRMRD::AcquisitionHeader acq_head = sptr_acq->getHead();
			
			if( acq_head.idx.contrast == i_contrast )
				acq_vec.append_acquisition_sptr( sptr_acq );
		}
		
		// std::shared_ptr< AcquisitionsVector > curr_template_acquis( new AcquisitionsVector(acq_vec) );
		auto curr_template_acquis = std::make_shared< AcquisitionsVector >(acq_vec);
		
		this->acq_model_.set_acquisition_template( curr_template_acquis );
		this->acq_model_.fwd(curr_img_wrap, csm_as_img, this->target_acquisitions_, offset);
		
	}
}

void PETDynamicSimulation::write_simulation_results( const std::string& filename_output_with_extension )
{
	
	std::cout << "Writing PET rawdata ... ";
	this->sptr_target_acquisitions_->write( filename_output_with_extension.c_str() );
	std::cout << "finished." << std::endl;

	this->pet_cont_gen_.map_tissue();

	std::vector< STIRImageData > contrast_filled_volumes = this->pet_cont_gen_.get_contrast_filled_volumes();

	STIRImageData attenuation_map = contrast_filled_volumes[1];

	attenuation_map = this->get_reduced_pet_img_in_template_format( attenuation_map );	

	std::stringstream stream_filename_attenuation_map; 
	stream_filename_attenuation_map << filename_output_with_extension.substr(0, filename_output_with_extension.length()-3);
	stream_filename_attenuation_map << "_attenuation_map.hv";

	attenuation_map.write( stream_filename_attenuation_map.str() );

	// stir::shared_ptr< stir::OutputFileFormat<Image3DF >> format_sptr =
	// stir::OutputFileFormat<Image3DF>::default_sptr();


	// Image3DF& attenuation_image = attenuation_map.data();

	// std::cout << "Writing PET attenuation map ... ";
	// format_sptr->write_to_file( stream_filename_attenuation_map.str() , attenuation_image);
	// std::cout << "... finished." << std::endl;
}

void PETDynamicSimulation::save_ground_truth_displacements( void )
{
	for(size_t i=0; i<this->motion_dynamics_.size(); i++)
	{
		this->motion_dynamics_[i]->save_ground_truth_displacements();
	}
}

void PETDynamicSimulation::add_dynamic( std::shared_ptr<PETMotionDynamic> sptr_motion_dyn)
{
	this->motion_dynamics_.push_back( sptr_motion_dyn );

}
void PETDynamicSimulation::add_dynamic( std::shared_ptr<PETContrastDynamic> sptr_contrast_dyn)
{
	this->contrast_dynamics_.push_back( sptr_contrast_dyn );
} 


void PETDynamicSimulation::simulate_statics()
{
	this->pet_cont_gen_.map_tissue();
	this->set_template_acquisition_data();
	this->acquire_raw_data();
	float const scale_factor = 25000.f;
	float const ms_per_second = 1000.f;
    float const result = scale_factor/ms_per_second;
    float const zero = 0.f;
	sptr_target_acquisitions_->axpby(&result, *sptr_target_acquisitions_, &zero, *sptr_target_acquisitions_ );

	this->add_noise();
	

}

void PETDynamicSimulation::add_noise( void )
{
	this->noise_generator_.add_noise( *sptr_target_acquisitions_, *sptr_target_acquisitions_ );
}	

void PETDynamicSimulation::simulate_dynamics( void )
{
	throw std::runtime_error( "Please give an acquisition time to simulate pet dynamics." );
}

void PETDynamicSimulation::simulate_dynamics( size_t const total_scan_time )
{
	if( motion_dynamics_.size() > 0)
		this->simulate_motion_dynamics( total_scan_time );		
	else
		throw std::runtime_error( "Please only use motion dynamics with PET for the moment." );
}

void PETDynamicSimulation::simulate_motion_dynamics(size_t const total_scan_time )
{
	cout << "Simulating PET motion dynamics... " <<endl;

	TimeBin total_time( 0, total_scan_time );
	TimeBinSet tot_time_interval{total_time};

	this->pet_cont_gen_.map_tissue();
	this->set_template_acquisition_data();

	size_t const num_motion_dynamics = this->motion_dynamics_.size();

	std::vector< int > all_num_dyn_states;
	for(size_t i=0; i<num_motion_dynamics; i++)
	{	
		cout << "Number simulated states for motion dynamic #" << i << ": " << motion_dynamics_[i]->get_num_simul_states() <<endl;
		all_num_dyn_states.push_back(motion_dynamics_[i]->get_num_simul_states());			
		motion_dynamics_[i]->prep_displacement_fields();
		motion_dynamics_[i]->align_motion_fields_with_image( this->template_image_data_);
	}

	LinearCombiGenerator lcg(all_num_dyn_states);
	
	size_t const num_total_dyn_states = lcg.get_num_total_combinations();
	std::vector< DimensionsType >  all_dyn_state_combos = lcg.get_all_combinations();
	
	

	for( size_t i_dyn_state=0; i_dyn_state < num_total_dyn_states; i_dyn_state++)
	{
		std::stringstream output_name_stream;
		output_name_stream << this->output_filename_prefix_;

		cout << "Acquisition motion state #" << i_dyn_state+1 << "/" << num_total_dyn_states << endl;

		DimensionsType current_combination = all_dyn_state_combos[i_dyn_state];

		TimeBinSet acquisition_time_bins_for_this_state = tot_time_interval;
		for( int i_motion_dyn = 0; i_motion_dyn<num_motion_dynamics; i_motion_dyn++ )
		{
			output_name_stream << "_dynamic_" << i_motion_dyn << "_state_" << current_combination[i_motion_dyn];	

			auto sptr_motion_dyn = this->motion_dynamics_[i_motion_dyn];
							
			TimeBinSet acquisition_times_in_bin = sptr_motion_dyn->get_time_bin_set_for_state( current_combination[i_motion_dyn] );
			acquisition_time_bins_for_this_state = intersect_time_bin_sets(acquisition_time_bins_for_this_state, acquisition_times_in_bin);

		}

		TimeAxisType time_in_dynamic_state = get_total_time_in_set(acquisition_time_bins_for_this_state);
		cout << "Time spent in this motion state: " << time_in_dynamic_state << endl;

		if( time_in_dynamic_state > 0)
		{
			std::vector<sirf::NiftiImageData3DDeformation<float> > all_motion_fields;
			
			for( int i_motion_dyn = 0; i_motion_dyn<num_motion_dynamics; i_motion_dyn++ )
			{
				cout << i_motion_dyn << endl;

				auto sptr_motion_dyn = this->motion_dynamics_[i_motion_dyn];
				std::vector< SignalBin > signal_bins = sptr_motion_dyn->get_bins();

				SignalBin bin = signal_bins[ current_combination[i_motion_dyn] ];	

				all_motion_fields.push_back( sptr_motion_dyn->get_interpolated_deformation_field( std::get<1>(bin) ) ); 

			}

			this->pet_cont_gen_.map_tissue();//crucial call, as the deformation results in deformed contrast generator data
			DynamicSimulationDeformer::deform_contrast_generator(this->pet_cont_gen_, all_motion_fields);
			this->acquire_raw_data();	

			float const ms_per_second = 1000.f;
            float const result = time_in_dynamic_state/ms_per_second;
            float const zero = 0.f;
			sptr_target_acquisitions_->axpby(&result, *sptr_target_acquisitions_, &zero, *sptr_target_acquisitions_ );

			this->add_noise();

			output_name_stream << ".hs";
			this->write_simulation_results( output_name_stream.str() );
		}
	}
}
		

void PETDynamicSimulation::set_template_acquisition_data(void)
{
	this->source_acquisitions_ = PETAcquisitionDataInFile( this->filename_rawdata_.c_str() );
}

void PETDynamicSimulation::set_template_image_data( const std::string& filename_header_with_ext )
{
	this->template_image_data_ = STIRImageData(filename_header_with_ext);
}

void PETDynamicSimulation::set_output_filename_prefix( const std::string& output_filename_prefix)
{
	this->output_filename_prefix_ = output_filename_prefix;
}

void PETDynamicSimulation::acquire_raw_data( void )
{
	std::vector< STIRImageData > contrast_filled_volumes = this->pet_cont_gen_.get_contrast_filled_volumes();

	STIRImageData activity_img = contrast_filled_volumes[0];
	STIRImageData attenuation_map = contrast_filled_volumes[1];

	activity_img = this->get_reduced_pet_img_in_template_format( activity_img );	
	attenuation_map = this->get_reduced_pet_img_in_template_format( attenuation_map );	


	STIRImageData template_img(activity_img);

	std::string outname_img = "/media/sf_SharedFolder/CCPPETMR/imageInMemory.hv";
	template_img.write( outname_img ); 

	stir::shared_ptr<stir::ProjMatrixByBin> sptr_ray_matrix (new RayTracingMatrix() );
	this->acq_model_.set_matrix( sptr_ray_matrix );		

	PETAttenuationModel att_mod(attenuation_map, this->acq_model_);
	this->acq_model_.set_asm( std::make_shared<PETAttenuationModel>(att_mod)); 

	auto succeeded = this->acq_model_.set_up( stir::shared_ptr<PETAcquisitionDataInFile>(new PETAcquisitionDataInFile(source_acquisitions_)),
	 			       stir::shared_ptr<STIRImageData>(new STIRImageData(template_img) ) );

	if( succeeded == stir::Succeeded::no )
		throw std::runtime_error("Setup of acquisition model failed");

	this->sptr_target_acquisitions_ = this->acq_model_.forward(activity_img);

	
}

STIRImageData PETDynamicSimulation::get_reduced_pet_img_in_template_format( const STIRImageData& full_size_img)
{
	std::vector< int > input_dims;
	input_dims.resize(3,0);
	full_size_img.get_dimensions(&input_dims[0]);

	size_t Nz = input_dims[2];
	size_t Ny = input_dims[1];
	size_t Nx = input_dims[0];
	
	

	size_t const num_voxels = Nx*Ny*Nz;

  	std::vector < float > vol_data;
  	vol_data.resize(num_voxels, 0);
  	full_size_img.get_data(&vol_data[0]);


	std::vector< int > template_dims;
	template_dims.resize(3,0);

	auto works = this->template_image_data_.get_dimensions(&template_dims[0]);

	if(works == -1)
		throw std::runtime_error("Irregular range of dimensions in PET image data.");

	std::reverse( input_dims.begin(), input_dims.end() );
	std::reverse( template_dims.begin(), template_dims.end() );

	std::vector< float > reduced_data;
	reduced_data.resize(template_dims[0]*template_dims[1]*template_dims[2],0);

	std::vector< size_t > offsets;
	for(int i = 0; i<3; i++)
	{
		if(input_dims[i] >= template_dims[i])
			offsets.push_back( size_t(float(input_dims[i]- template_dims[i])/2.f));
	
		else
			throw std::runtime_error("Please give only data which has equal or larger data dimensions than the template image.");
	}

	offsets[2] = input_dims[2] - template_dims[2];

	for(size_t nz = 0; nz<template_dims[2]; nz++)
	for(size_t ny = 0; ny<template_dims[1]; ny++)
	for(size_t nx = 0; nx<template_dims[0]; nx++)
	{
		
		size_t const linear_index_vol_data = ( (nz+offsets[2]) * input_dims[1] + (ny+offsets[1]) ) * input_dims[0] + (nx+offsets[0]);
		size_t const linear_index_subset = (nz*template_dims[1] + ny)*template_dims[0] + nx;
		
		reduced_data[linear_index_subset] = vol_data[linear_index_vol_data];
	}	

	STIRImageData out( this-> template_image_data_ );
	out.set_data(&reduced_data[0]);

	return out;

}