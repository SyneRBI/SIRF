/* ================================================

Author: Johannes Mayer
Date: 2018.07.20
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */
#pragma once

#include <string>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "gadgetron_data_containers.h"


#include "gadgetron_x.h"
#include "tissuelabelmapper.h"
#include "contrastgenerator.h"



class MRDynamicSimulation {


public:

	MRDynamicSimulation( MRContrastGenerator mr_cont_gen) : mr_cont_gen_( mr_cont_gen ) { };
	
	void set_filename_mr_rawdata( std::string const filename_ismrmrd_src_file ) { this->filename_mr_rawdata_ = filename_ismrmrd_src_file; };
	std::string get_filename_mr_rawdata( void )	{ return this-> filename_mr_rawdata_; };

	ISMRMRD::IsmrmrdHeader get_ismrmrd_header( void ){ return this->hdr_;};

	void simulate_dynamics( void );

	void write_simulation_results( std::string const filename_output_with_h5_extension );

	void extract_src_information( void );

protected:
	
	std::string filename_mr_rawdata_;
	ISMRMRD::IsmrmrdHeader hdr_;

	AcquisitionsVector source_acquisitions_;
	AcquisitionsVector target_acquisitions_;
	
	MRContrastGenerator mr_cont_gen_;
	//MRAcquisitionModel acq_model_;


};


