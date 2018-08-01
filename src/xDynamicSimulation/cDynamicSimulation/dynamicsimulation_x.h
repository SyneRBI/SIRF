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


#define IMG_DATA_TYPE 7 // from ismrmrd enum ISMRMRD_CXFLOAT = 7





class aDynamicSimulation {

public:

	aDynamicSimulation(){};
	~aDynamicSimulation(){};

	virtual std::string get_filename_rawdata( void )
	{ 
		return this-> filename_rawdata_; 
	}

	virtual	void set_filename_rawdata( std::string const filename_template_rawdata ) 
	{ 
		this->filename_rawdata_ = filename_template_rawdata; 
	}

	virtual void simulate_dynamics( void ) = 0;
	virtual void write_simulation_results( std::string const filename_output_with_extension ) = 0;
	
	/*void add_dynamic( MotionDynamic motion_dyn) 
	{
		this->motion_dynamics_.push_back(motion_dyn);
	};
	void add_dynamic( ContrastDynamic cont_dyn) 
	{
		this->contrast_dynamics_.push_back(cont_dyn);
	};*/


protected:

	/*std::vector< MotionDynamic > motion_dynamics_;
	std::vector< ContrastDynamic > contrast_dynamics_;*/

	std::string filename_rawdata_;

};




class MRDynamicSimulation : public aDynamicSimulation {

public:

	MRDynamicSimulation( MRContrastGenerator mr_cont_gen) : mr_cont_gen_( mr_cont_gen ) { };
	void write_simulation_results( std::string const filename_output_with_extension );


	ISMRMRD::IsmrmrdHeader get_ismrmrd_header( void ){ return this->hdr_;};
	void extract_src_information( void );

	void simulate_dynamics( void );

private:

	ISMRMRD::IsmrmrdHeader hdr_;

	AcquisitionsVector source_acquisitions_;
	AcquisitionsVector target_acquisitions_;
	
	MRContrastGenerator mr_cont_gen_;
	MRAcquisitionModel acq_model_;

};


/*class PETDynamicSimulation : public aDynamicSimulation{

public:
	PETDynamicSimulation( PETContrastGenerator pet_cont_gen) : pet_cont_gen_(pet_cont_gen);

	void simulate_dynamics( void ){};

private:

	PETContrastGenerator pet_cont_gen_;
	
};
*/

