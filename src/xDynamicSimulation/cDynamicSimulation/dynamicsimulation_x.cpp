/*
author: johannes mayer
date: 15. March 2018

*/


#include "dynamicsimulation_x.h"

#include "auxiliary_input_output.h"


void MRDynamicSimulation::write_simulation_results( std::string const filename_output_with_h5_extension ) 
{
	target_acquisitions_.write( filename_output_with_h5_extension.c_str() );
}

void MRDynamicSimulation::simulate_dynamics( void )
{

	this->extract_src_information();
	this->mr_cont_gen_.map_contrast();
	std::vector< ISMRMRD::Image< complex_float_t> > contrast_filled_volumes = this->mr_cont_gen_.get_contrast_filled_volumes();

	size_t const num_contrasts = contrast_filled_volumes.size();

	std::cout << num_contrasts << std::endl;

}

void MRDynamicSimulation::extract_src_information( void )
{
	this->hdr_ = mr_io::read_ismrmrd_header( filename_mr_rawdata_ );

	this->mr_cont_gen_.set_rawdata_header( this->hdr_ );

	this->source_acquisitions_ = mr_io::read_ismrmrd_acquisitions( filename_mr_rawdata_ );
	this->target_acquisitions_.copy_acquisitions_info( this->source_acquisitions_ );

}



