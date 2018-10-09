/*
author: johannes mayer
date: 15. March 2018

*/

#include <memory>
#include <mutex>
#include <sstream>

#include "dynamicsimulation_x.h"

#include "auxiliary_input_output.h"

#include "SIRFImageDataDeformation.h"

#include "dynsim_deformer.h"



using namespace sirf;
using std::cout;
using std::endl;

void MRDynamicSimulation::write_simulation_results( std::string const filename_output_with_h5_extension ) 
{	
	try	
	{
		cout << "Started writing simulation output to: " << filename_output_with_h5_extension <<endl;
		cout << "Number of acquisitions to write: " << this->target_acquisitions_.number() << endl;

		std::stringstream serialized_hdr;
		ISMRMRD::serialize(this->hdr_, serialized_hdr);
		target_acquisitions_.set_acquisitions_info( serialized_hdr.str() ); 

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

void MRDynamicSimulation::add_dynamic( std::shared_ptr<MRMotionDynamic> sptr_motion_dyn )
{
	this->motion_dynamics_.push_back( sptr_motion_dyn );
};


void MRDynamicSimulation::add_dynamic( std::shared_ptr<MRContrastDynamic> sptr_contrast_dyn) 
{
	this->contrast_dynamics_.push_back( sptr_contrast_dyn );
};


void MRDynamicSimulation::simulate_statics( void )
{
	cout << "Simulating static data acquisition... " <<endl;



	this->extract_hdr_information();

	this->acq_model_.setTraj( this->sptr_trajectory_ );

	this->mr_cont_gen_.map_contrast();
	this->source_acquisitions_ = this->all_source_acquisitions_;
	this->acquire_raw_data();
}

void MRDynamicSimulation::simulate_dynamics( void )
{
	if( this->motion_dynamics_.size() > 0 && this->contrast_dynamics_.size() > 0 )
		this->simulate_simultaneous_motion_contrast_dynamics();		
	else if( this->motion_dynamics_.size() > 0 && this->contrast_dynamics_.size() == 0 )
		this->simulate_motion_dynamics();
	else if (this->motion_dynamics_.size() == 0 && this->contrast_dynamics_.size() > 0)
		this->simulate_contrast_dynamics();
	else
		cout << "No dynamics added to the simulation. Well, if you don't give any dynamics, simulating them is not possible." << endl;
}

void MRDynamicSimulation::simulate_simultaneous_motion_contrast_dynamics()
{
	cout << "Simulating motion and contrast dynamics... " <<endl;

	this->extract_hdr_information();
	this->acq_model_.setTraj( this->sptr_trajectory_ );

	this->mr_cont_gen_.map_contrast();
		
	size_t const num_contrast_dyns = this->contrast_dynamics_.size();
	size_t const num_motion_dyns = this->motion_dynamics_.size();

	std::vector< int > num_all_dyn_states;

	for(size_t i=0; i<num_contrast_dyns; i++)
		num_all_dyn_states.push_back(contrast_dynamics_[i]->get_num_simul_states());			

	for(size_t i=0; i<num_motion_dyns; i++)	
	{
		num_all_dyn_states.push_back(motion_dynamics_[i]->get_num_simul_states());			
		motion_dynamics_[i]->prep_displacements_fields();
	}

	LinearCombiGenerator lcg(num_all_dyn_states);
	
	size_t const num_total_dyn_states = lcg.get_num_total_combinations();

	std::vector< DimensionsType >  all_dynamic_state_combos = lcg.get_all_combinations();

	for( size_t i_dyn_state=0; i_dyn_state < num_total_dyn_states; i_dyn_state++)
	{
		cout << "Acquisition dynamic state #" << i_dyn_state << "/" << num_total_dyn_states << endl;

		DimensionsType current_combination = all_dynamic_state_combos[i_dyn_state];

		sirf::AcquisitionsVector acquisitions_for_this_state = this->all_source_acquisitions_;

		
		for( int i_contrast_dyn = 0; i_contrast_dyn<num_contrast_dyns; i_contrast_dyn++ )
		{
			auto sptr_contrast_dyn = this->contrast_dynamics_[i_contrast_dyn];


			int const current_bin = current_combination[i_contrast_dyn];						

			AcquisitionsVector acquis_in_bin = sptr_contrast_dyn->get_binned_mr_acquisitions( current_bin );
			acquisitions_for_this_state = intersect_mr_acquisition_data(acquisitions_for_this_state, acquis_in_bin);

		}

		for( int i_motion_dyn = 0; i_motion_dyn<num_motion_dyns; i_motion_dyn++ )
		{
			auto motion_dyn = this->motion_dynamics_[ i_motion_dyn ];
			int const current_bin = current_combination[ num_contrast_dyns + i_motion_dyn ];						

			AcquisitionsVector acquis_in_bin = motion_dyn->get_binned_mr_acquisitions( current_bin );
			acquisitions_for_this_state = intersect_mr_acquisition_data(acquisitions_for_this_state, acquis_in_bin);

		}

		cout << "# of mr acquis in this dynamic state: " << acquisitions_for_this_state.number() << endl;
		
		if( acquisitions_for_this_state.number() > 0)
		{
			for( int i_contrast_dyn = 0; i_contrast_dyn<num_contrast_dyns; i_contrast_dyn++ )
			{
				auto sptr_contrast_dyn = this->contrast_dynamics_[i_contrast_dyn];

				std::vector< SignalBin > signal_bins = sptr_contrast_dyn->get_bins();

				int const contrast_bin_number = current_combination[i_contrast_dyn];						
				SignalBin bin = signal_bins[ contrast_bin_number ];	

				TissueParameterList tissueparameter_list_to_replace = sptr_contrast_dyn->get_interpolated_tissue_params( std::get<1>(bin) );

				for( size_t i_tiss=0; i_tiss< tissueparameter_list_to_replace.size(); i_tiss++ )
				{
					TissueParameter curr_param = tissueparameter_list_to_replace[i_tiss];
					this->mr_cont_gen_.replace_petmr_tissue_parameters( curr_param.label_, curr_param );	
				}
			}

			std::vector<SIRFImageDataDeformation> all_motion_fields;
			
			for( int i_motion_dyn = 0; i_motion_dyn<num_motion_dyns; i_motion_dyn++ )
			{
				cout << i_motion_dyn << endl;

				auto motion_dyn = this->motion_dynamics_[i_motion_dyn];
				std::vector< SignalBin > signal_bins = motion_dyn->get_bins();

				int const motion_bin_number = current_combination[ num_contrast_dyns + i_motion_dyn ];						
				SignalBin bin = signal_bins[ motion_bin_number ];

				all_motion_fields.push_back( motion_dyn->get_interpolated_displacement_field( std::get<1>(bin) ) ); 

			}

			this->mr_cont_gen_.map_contrast();//crucial call, as the deformation results in deformed contrast generator data
			cout << " num motoin fields " << all_motion_fields.size() << endl;
			DynamicSimulationDeformer::deform_contrast_generator(this->mr_cont_gen_, all_motion_fields);
			this->source_acquisitions_ = acquisitions_for_this_state;
			this->acquire_raw_data();	
		}

	}
	this->noise_generator_.add_noise(this->target_acquisitions_);

	// for(size_t i=0; i<num_motion_dyns; i++)
	// 	this->motion_dynamics_[i].delete_temp_folder();	

}
void MRDynamicSimulation::simulate_contrast_dynamics( void )
{
	cout << "Simulating contrast dynamics... " <<endl;

	this->extract_hdr_information();
	this->acq_model_.setTraj( this->sptr_trajectory_ );

	this->mr_cont_gen_.map_contrast();

		
	size_t const num_contrast_dyns = this->contrast_dynamics_.size();

	std::vector< int > all_num_dyn_states;
	for(size_t i=0; i<num_contrast_dyns; i++)
		all_num_dyn_states.push_back(contrast_dynamics_[i]->get_num_simul_states());			


	LinearCombiGenerator lcg(all_num_dyn_states);
	
	size_t const num_total_dyn_states = lcg.get_num_total_combinations();
	std::vector< DimensionsType >  all_dyn_state_combos = lcg.get_all_combinations();
	
	for( size_t i_dyn_state=0; i_dyn_state < num_total_dyn_states; i_dyn_state++)
	{
		cout << "Acquisition dynamic state #" << i_dyn_state << "/" << num_total_dyn_states << endl;

		DimensionsType current_combination = all_dyn_state_combos[i_dyn_state];

		sirf::AcquisitionsVector acquisitions_for_this_state = this->all_source_acquisitions_;


		for( int i_contrast_dyn = 0; i_contrast_dyn<num_contrast_dyns; i_contrast_dyn++ )
		{
			auto sptr_contrast_dyn = this->contrast_dynamics_[i_contrast_dyn];
			std::vector< SignalBin > signal_bins = sptr_contrast_dyn->get_bins();

			SignalBin bin = signal_bins[ current_combination[i_contrast_dyn] ];	
			TissueParameterList tissueparameter_list_to_replace = sptr_contrast_dyn->get_interpolated_tissue_params( std::get<1>(bin) );

			for( size_t i_tiss=0; i_tiss< tissueparameter_list_to_replace.size(); i_tiss++ )
			{
				TissueParameter curr_param = tissueparameter_list_to_replace[i_tiss];
				this->mr_cont_gen_.replace_petmr_tissue_parameters( curr_param.label_, curr_param );	
			}

			AcquisitionsVector acquis_in_bin = sptr_contrast_dyn->get_binned_mr_acquisitions( current_combination[i_contrast_dyn] );
			acquisitions_for_this_state = intersect_mr_acquisition_data(acquisitions_for_this_state, acquis_in_bin);

		}

		cout << "# of acquis in this dynamic state: " << acquisitions_for_this_state.number() << endl;

		if( acquisitions_for_this_state.number() > 0)
		{
			this->mr_cont_gen_.map_contrast();
			this->source_acquisitions_ = acquisitions_for_this_state;
			this->acquire_raw_data();	
		}
	}

	this->noise_generator_.add_noise(this->target_acquisitions_);


}

void MRDynamicSimulation::simulate_motion_dynamics( void )
{
	cout << "Simulating motion dynamics... " <<endl;

	this->extract_hdr_information();
	this->acq_model_.setTraj( this->sptr_trajectory_ );

	this->mr_cont_gen_.map_contrast();

	size_t const num_motion_dynamics = this->motion_dynamics_.size();

	std::vector< int > all_num_dyn_states;
	for(size_t i=0; i<num_motion_dynamics; i++)
	{	
		cout << motion_dynamics_[i]->get_num_simul_states() <<endl;
		all_num_dyn_states.push_back(motion_dynamics_[i]->get_num_simul_states());			
		motion_dynamics_[i]->prep_displacements_fields();
	}

	LinearCombiGenerator lcg(all_num_dyn_states);
	
	size_t const num_total_dyn_states = lcg.get_num_total_combinations();
	std::vector< DimensionsType >  all_dyn_state_combos = lcg.get_all_combinations();
	

	for( size_t i_dyn_state=0; i_dyn_state < num_total_dyn_states; i_dyn_state++)
	{
		cout << "Acquisition motion state #" << i_dyn_state << "/" << num_total_dyn_states << endl;

		DimensionsType current_combination = all_dyn_state_combos[i_dyn_state];

		sirf::AcquisitionsVector acquisitions_for_this_state = this->all_source_acquisitions_;

		for( int i_motion_dyn = 0; i_motion_dyn<num_motion_dynamics; i_motion_dyn++ )
		{
			auto sptr_motion_dyn = this->motion_dynamics_[i_motion_dyn];
						
			AcquisitionsVector acquis_in_bin = sptr_motion_dyn->get_binned_mr_acquisitions( current_combination[i_motion_dyn] );
			acquisitions_for_this_state = intersect_mr_acquisition_data(acquisitions_for_this_state, acquis_in_bin);

		}

		cout << "# of mr acquis in this motion state: " << acquisitions_for_this_state.number() << endl;

		if( acquisitions_for_this_state.number() > 0)
		{
			std::vector<SIRFImageDataDeformation> all_motion_fields;
			
			for( int i_motion_dyn = 0; i_motion_dyn<num_motion_dynamics; i_motion_dyn++ )
			{
				cout << i_motion_dyn << endl;

				auto sptr_motion_dyn = this->motion_dynamics_[i_motion_dyn];
				std::vector< SignalBin > signal_bins = sptr_motion_dyn->get_bins();

				SignalBin bin = signal_bins[ current_combination[i_motion_dyn] ];	

				all_motion_fields.push_back( sptr_motion_dyn->get_interpolated_displacement_field( std::get<1>(bin) ) ); 

			}

			this->mr_cont_gen_.map_contrast();//crucial call, as the deformation results in deformed contrast generator data
			DynamicSimulationDeformer::deform_contrast_generator(this->mr_cont_gen_, all_motion_fields);
			this->source_acquisitions_ = acquisitions_for_this_state;
			this->acquire_raw_data();	
		}
	}

	this->noise_generator_.add_noise(this->target_acquisitions_);

	for(size_t i=0; i<num_motion_dynamics; i++)
		motion_dynamics_[i]->delete_temp_folder();		
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


void MRDynamicSimulation::set_trajectory( std::shared_ptr<sirf::aTrajectoryContainer> sptr_trajectory)
{
	this->sptr_trajectory_ = sptr_trajectory;
	sptr_trajectory_->overwrite_ismrmrd_trajectory_info( this->hdr_ );
	
	this->mr_cont_gen_.set_rawdata_header( this->hdr_ );

}


void MRDynamicSimulation::set_all_source_acquisitions(MRDataContainerType acquisitions )
{
	this->all_source_acquisitions_ = acquisitions;
	this->target_acquisitions_.copy_acquisitions_info( this->all_source_acquisitions_ );
}


void MRDynamicSimulation::set_noise_width(float const sigma)
{
	this->noise_generator_.set_noise_width( sigma );
}

void MRDynamicSimulation::set_SNR(float const SNR)
{
	this->noise_generator_.set_SNR(SNR);
}

void MRDynamicSimulation::acquire_raw_data( void )
{

	std::vector< ISMRMRD::Image< complex_float_t> > contrast_filled_volumes = this->mr_cont_gen_.get_contrast_filled_volumes();

	size_t const num_contrasts = contrast_filled_volumes.size();

	size_t Nx = contrast_filled_volumes[0].getMatrixSizeX();
	size_t Ny = contrast_filled_volumes[0].getMatrixSizeY();
	size_t Nz = contrast_filled_volumes[0].getMatrixSizeZ();

	size_t Nc = 4;
	CoilDataAsCFImage csm_as_img( Nx, Ny, Nz , Nc);
	std::vector< complex_float_t > mock_csm;
	mock_csm.resize( Nx * Ny * Nz * Nc, std::complex<float>(1,0) );

	csm_as_img.set_data( &mock_csm[0] );

	unsigned int offset = 0;


	for( size_t i_contrast=0; i_contrast<num_contrasts; i_contrast++)
	{
		ISMRMRD::Acquisition acq;

		cout << "Acquisition contrast " << i_contrast << endl;
		ISMRMRD::Image<complex_float_t> curr_cont = contrast_filled_volumes[i_contrast];
		curr_cont = vol_orientator_.reorient_image(curr_cont);

		ImageWrap curr_img_wrap(IMG_DATA_TYPE, new ISMRMRD::Image< complex_float_t >(curr_cont));		

		AcquisitionsVector acq_vec;
		acq_vec.copy_acquisitions_info( this->source_acquisitions_ );
		
		for( size_t i_acqu=0; i_acqu<this->source_acquisitions_.items(); i_acqu++)
		{
			this->source_acquisitions_.get_acquisition(i_acqu, acq);
			ISMRMRD::AcquisitionHeader acq_head = acq.getHead();
			
			if( acq_head.idx.contrast == i_contrast )
				acq_vec.append_acquisition( acq );

		}
		
		std::shared_ptr< AcquisitionsVector > curr_template_acquis( new AcquisitionsVector(acq_vec) );
		
		this->acq_model_.set_acquisition_template( curr_template_acquis );
		
		this->acq_model_.fwd(curr_img_wrap, csm_as_img, this->target_acquisitions_, offset);
		
	}
}

void PETDynamicSimulation::write_simulation_results( std::string const filename_output_with_extension )
{
	this->target_acquisitions_->write( filename_output_with_extension.c_str() );
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

	target_acquisitions_
}
		

void PETDynamicSimulation::simulate_dynamics()
{
	this->pet_cont_gen_.map_tissue();
	this->set_template_acquisition_data();
	this->acquire_raw_data();
}
		

void PETDynamicSimulation::set_template_acquisition_data(void)
{
	this->source_acquisitions_ = PETAcquisitionDataInFile( this->filename_rawdata_.c_str() );
}

void PETDynamicSimulation::set_template_image_data( std::string const filename_header_with_ext )
{
	this->template_image_data_ = sirf::PETImageData(filename_header_with_ext);
}


void PETDynamicSimulation::acquire_raw_data( void )
{
	std::vector< sirf::PETImageData > contrast_filled_volumes = this->pet_cont_gen_.get_contrast_filled_volumes();

	sirf::PETImageData activity_img = contrast_filled_volumes[0];
	sirf::PETImageData attenaution_map = contrast_filled_volumes[1];

	activity_img = this->get_reduced_pet_img_in_template_format( activity_img );	
	attenaution_map = this->get_reduced_pet_img_in_template_format( attenaution_map );	

	sirf::PETImageData template_img(activity_img);

	sirf::Image3DF& image = activity_img.data();
	stir::shared_ptr< stir::OutputFileFormat<sirf::Image3DF >> format_sptr =
	stir::OutputFileFormat<sirf::Image3DF>::default_sptr();

	format_sptr->write_to_file( "/media/sf_SharedFolder/CCPPETMR/imageInMemory.hv" , image);
	cout << "Finished Writing image in memory" << endl;
	
	stir::shared_ptr<stir::ProjMatrixByBin> sptr_ray_matrix (new sirf::RayTracingMatrix() );

	this->acq_model_.set_matrix( sptr_ray_matrix );		

	auto succeeded = this->acq_model_.set_up( stir::shared_ptr<PETAcquisitionDataInFile>(new PETAcquisitionDataInFile(source_acquisitions_)),
	 			       stir::shared_ptr<PETImageData>(new PETImageData(template_img) ) );	

	if( succeeded == stir::Succeeded::no )
		throw std::runtime_error("Setup of acquisition model failed");

		
	this->target_acquisitions_ = this->acq_model_.forward(activity_img);

	
}

sirf::PETImageData PETDynamicSimulation::get_reduced_pet_img_in_template_format( const sirf::PETImageData& full_size_img)
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

	std::reverse( template_dims.begin(), template_dims.end() );

	std::vector< float > reduced_data;
	reduced_data.resize(template_dims[0]*template_dims[1]*template_dims[2],0);

	std::vector< size_t > offsets;
	for(int i = 0; i<3; i++)
	{
		if(input_dims[i] >= template_dims[i])
			offsets.push_back(input_dims[i]/2 - template_dims[i]/2);
		else
			throw std::runtime_error("Please give only data which has equal or larger data dimensions than the template image.");
	}

	for(size_t nz = 0; nz<template_dims[2]; nz++)
	for(size_t ny = 0; ny<template_dims[1]; ny++)
	for(size_t nx = 0; nx<template_dims[0]; nx++)
	{
		
		size_t const linear_index_vol_data = ( (nz+offsets[2]) * input_dims[1] + (ny+offsets[1]) ) * input_dims[0] + (nx+offsets[0]);
		size_t const linear_index_subset = (nz*template_dims[1] + ny)*template_dims[0] + nx;
		
		reduced_data[linear_index_subset] = vol_data[linear_index_vol_data];
	}	

	PETImageData out( this-> template_image_data_ );
	out.set_data(&reduced_data[0]);

	return out;

}