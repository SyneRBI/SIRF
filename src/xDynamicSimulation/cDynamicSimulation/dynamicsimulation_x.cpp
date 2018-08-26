/*
author: johannes mayer
date: 15. March 2018

*/


#include "dynamicsimulation_x.h"

#include "auxiliary_input_output.h"

#include <memory>

using namespace sirf;

void MRDynamicSimulation::write_simulation_results( std::string const filename_output_with_h5_extension ) 
{	
	try	
	{
		std::cout << "Started writing simulation output to: " << filename_output_with_h5_extension <<std::endl;
		std::cout << "Number of acquisitions to write: " << this->target_acquisitions_.number() << std::endl;
		target_acquisitions_.write( filename_output_with_h5_extension.c_str() );
		std::cout << "Finished writing simulation output."<<std::endl;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "Maybe you forgot to give the correct filename ending .h5 in this case." << std::endl;
	}
}

void MRDynamicSimulation::simulate_dynamics( void )
{
	this->extract_hdr_information();
	this->mr_cont_gen_.map_contrast();

	if(this->motion_dynamics_.size() != 0)
		throw std::runtime_error("So far only one contrast and zero motion dynamics are supported. Please give the appropriate number of dynamics.");			

	
	size_t const num_contrast_dyns = this->contrast_dynamics_.size();

	std::vector< int > all_num_dyn_states;
	for(size_t i=0; i<num_contrast_dyns; i++)
		all_num_dyn_states.push_back(contrast_dynamics_[i].get_num_simul_states());			


	LinearCombiGenerator lcg(all_num_dyn_states);
	
	size_t const num_total_dyn_states = lcg.get_num_total_combinations();
	std::vector< DimensionsType >  all_dyn_state_combos = lcg.get_all_combinations();
	
	for( size_t i_dyn_state=0; i_dyn_state < num_total_dyn_states; i_dyn_state++)
	{
		std::cout << "Acquiriing dynamic state #" << i_dyn_state << "/" << num_total_dyn_states << std::endl;

		DimensionsType current_combination = all_dyn_state_combos[i_dyn_state];

		sirf::AcquisitionsVector acquisitions_for_this_state = this->all_source_acquisitions_;

		std::vector< SignalAxisType> cont_signals;

		for( int i_cont_dyn = 0; i_cont_dyn<num_contrast_dyns; i_cont_dyn++ )
		{
			ContrastDynamic cont_dyn = this->contrast_dynamics_[i_cont_dyn];
			std::vector< SignalBin > signal_bins = cont_dyn.get_bins();

			SignalBin bin = signal_bins[ current_combination[i_cont_dyn] ];	
			TissueParameterList tissueparameter_list_to_replace = cont_dyn.get_interpolated_tissue_params( std::get<1>(bin) );
			cont_signals.push_back(std::get<1>(bin));	

			for( size_t i_tiss=0; i_tiss< tissueparameter_list_to_replace.size(); i_tiss++ )
			{
				TissueParameter curr_param = tissueparameter_list_to_replace[i_tiss];
				this->mr_cont_gen_.replace_petmr_tissue_parameters( curr_param.label_, curr_param );	
			}

			AcquisitionsVector acquis_in_bin = cont_dyn.get_binned_mr_acquisitions( current_combination[i_cont_dyn] );
			acquisitions_for_this_state = intersect_mr_acquisition_data(acquisitions_for_this_state, acquis_in_bin);

		}

		std::cout << "# of acquis in this dynamic state: " << acquisitions_for_this_state.number() << std::endl;

		if( acquisitions_for_this_state.number() > 0)
		{
			this->mr_cont_gen_.map_contrast();
			this->source_acquisitions_ = acquisitions_for_this_state;
			this->acquire_raw_data();	
		}
	}

	this->noise_generator_.add_noise(this->target_acquisitions_);
}

void MRDynamicSimulation::extract_hdr_information( void )
{
	this->hdr_ = mr_io::read_ismrmrd_header( filename_rawdata_ );

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

		std::cout << "Acquisition contrast " << i_contrast << std::endl;
		ISMRMRD::Image<complex_float_t> curr_cont = contrast_filled_volumes[i_contrast];
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




void PETDynamicSimulation::simulate_dynamics()
{
	this->pet_cont_gen_.map_contrast();
	this->set_template_acquisition_data();
	this->acquire_raw_data();
}
		
void PETDynamicSimulation::set_template_acquisition_data(void)
{
	this->source_acquisitions_ = PETAcquisitionDataInFile( this->filename_rawdata_.c_str() );
}


void PETDynamicSimulation::acquire_raw_data( void )
{
	std::vector< sirf::PETImageData > contrast_filled_volumes = this->pet_cont_gen_.get_contrast_filled_volumes();
	size_t num_contrast = contrast_filled_volumes.size();

	if( num_contrast != 1)
		throw std::runtime_error("Please give only one contrast for the time being."); // POSSIBLY REMOVE LATER!

	auto dims = pet_cont_gen_.get_dimensions();
	auto vxsizes = pet_cont_gen_.get_voxel_sizes();

	for( size_t i_vol=0; i_vol<num_contrast; i_vol++)
	{

		sirf::PETImageData curr_img = contrast_filled_volumes[i_vol];
		sirf::PETImageData template_img(curr_img);

		// sirf::Image3DF& image = curr_img.data();
		// stir::shared_ptr< stir::OutputFileFormat<sirf::Image3DF >> format_sptr =
		// stir::OutputFileFormat<sirf::Image3DF>::default_sptr();

		// format_sptr->write_to_file( "/media/sf_SharedFolder/CCPPETMR/imageInMemory.hv" , image);
		// std::cout << "Finished Writing image in memory" << std::endl;
		
		stir::shared_ptr<stir::ProjMatrixByBin> sptr_ray_matrix (new sirf::RayTracingMatrix() );

		this->acq_model_.set_matrix( sptr_ray_matrix );		

		auto succeeded = this->acq_model_.set_up( stir::shared_ptr<PETAcquisitionDataInFile>(new PETAcquisitionDataInFile(source_acquisitions_)),
		 			       stir::shared_ptr<PETImageData>(new PETImageData(template_img) ) );	

		if( succeeded == stir::Succeeded::no )
			throw std::runtime_error("Setup of acquisition model failed");

			
		this->target_acquisitions_ = this->acq_model_.forward(curr_img);

	 }
}