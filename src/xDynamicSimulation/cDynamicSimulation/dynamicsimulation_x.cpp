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

	this->extract_src_information();
	this->mr_cont_gen_.map_contrast();
	//this->mr_cont_gen_.match_output_dims_to_headerinfo();


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

	ISMRMRD::Acquisition acq;
	
	for( size_t i_contrast=0; i_contrast<num_contrasts; i_contrast++)
	{
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

void MRDynamicSimulation::extract_src_information( void )
{
	this->hdr_ = mr_io::read_ismrmrd_header( filename_rawdata_ );

	this->mr_cont_gen_.set_rawdata_header( this->hdr_ );

	this->source_acquisitions_ = mr_io::read_ismrmrd_acquisitions( filename_rawdata_ );
	this->target_acquisitions_.copy_acquisitions_info( this->source_acquisitions_ );

}


void PETDynamicSimulation::write_simulation_results( std::string const filename_output_with_extension )
{
	this->target_acquisitions_->write( filename_output_with_extension.c_str() );
}


void PETDynamicSimulation::set_template_acquisition_data(void)
{
	this->source_acquisitions_ = PETAcquisitionDataInFile( this->filename_rawdata_.c_str() );
}

void PETDynamicSimulation::simulate_dynamics()
{

	this->pet_cont_gen_.map_contrast();
	std::vector< sirf::PETImageData > contrast_filled_volumes = this->pet_cont_gen_.get_contrast_filled_volumes();
	size_t num_contrast = contrast_filled_volumes.size();

	if( num_contrast != 1)
		throw std::runtime_error("Please give only one contrast for the time being."); // POSSIBLY REMOVE LATER!



	auto dims = pet_cont_gen_.get_dimensions();
	auto vxsizes = pet_cont_gen_.get_voxel_sizes();

/*
	std::cout << dims[0] << std::endl;
	std::cout << dims[1] << std::endl;
	std::cout << dims[2] << std::endl;

	std::cout << vxsizes[0] << std::endl;
	std::cout << vxsizes[1] << std::endl;
	std::cout << vxsizes[2] << std::endl;
*/
	for( size_t i_vol=0; i_vol<num_contrast; i_vol++)
	{

		sirf::PETImageData curr_img = contrast_filled_volumes[i_vol];
		sirf::PETImageData template_img(curr_img);

		sirf::Image3DF& image = curr_img.data();
		stir::shared_ptr< stir::OutputFileFormat<sirf::Image3DF >> format_sptr =
		stir::OutputFileFormat<sirf::Image3DF>::default_sptr();

		format_sptr->write_to_file( "/media/sf_SharedFolder/CCPPETMR/imageInMemory.hv" , image);

		
		stir::shared_ptr<stir::ProjMatrixByBin> sptr_ray_matrix (new sirf::RayTracingMatrix() );
		std::cout << "a" << std::endl;
		this->acq_model_.set_matrix( sptr_ray_matrix );		
		std::cout << "b" << std::endl;
		auto succeeded = this->acq_model_.set_up( stir::shared_ptr<PETAcquisitionDataInFile>(new PETAcquisitionDataInFile(source_acquisitions_)),
		 			       stir::shared_ptr<PETImageData>(new PETImageData(template_img) ) );	

		if( succeeded == stir::Succeeded::no )
			throw std::runtime_error("Setup of acquisition model failed");

			
		this->target_acquisitions_ = this->acq_model_.forward(curr_img);

	 }

}
		