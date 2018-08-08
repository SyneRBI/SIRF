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
#include "stir_x.h"

#include "tissuelabelmapper.h"
#include "contrastgenerator.h"
#include "dynamics.h"



#define IMG_DATA_TYPE 7 // from ismrmrd enum ISMRMRD_CXFLOAT = 7





class aDynamicSimulation {

public:

	aDynamicSimulation(){};
	~aDynamicSimulation(){};

	std::string get_filename_rawdata( void )
	{ 
		return this-> filename_rawdata_; 
	}

	void set_filename_rawdata( std::string const filename_template_rawdata ) 
	{ 
		this->filename_rawdata_ = filename_template_rawdata; 
	}

	virtual void simulate_dynamics( void ) = 0;
	virtual void write_simulation_results( std::string const filename_output_with_extension ) = 0;

	void add_dynamic( MotionDynamic motion_dyn) 
	{
		this->motion_dynamics_.push_back(motion_dyn);
	};
	/*
	void add_dynamic( ContrastDynamic cont_dyn) 
	{
		this->contrast_dynamics_.push_back(cont_dyn);
	};*/


protected:

	std::vector< MotionDynamic > motion_dynamics_;
	// std::vector< ContrastDynamic > contrast_dynamics_;

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

	sirf::AcquisitionsVector source_acquisitions_;
	sirf::AcquisitionsVector target_acquisitions_;
	
	MRContrastGenerator mr_cont_gen_;
	sirf::MRAcquisitionModel acq_model_;

};


class PETDynamicSimulation : public aDynamicSimulation{

public:
	PETDynamicSimulation( PETContrastGenerator pet_cont_gen ) : pet_cont_gen_(pet_cont_gen) 
	{ 
		// this->pet_cont_gen_ = pet_cont_gen;
	};
		
	void simulate_dynamics( void ){};
	void write_simulation_results( std::string const filename_output_with_extension );

	void extract_src_information( void );

private:

	PETContrastGenerator pet_cont_gen_;
	sirf::PETAcquisitionModelUsingMatrix acq_model_;

	sirf::PETAcquisitionDataInFile source_acquisitions_;
	sirf::PETAcquisitionDataInMemory target_acquisitions_;
	
};








class LinearCombiGenerator{

public:
	LinearCombiGenerator(std::vector<size_t> const dims): dims_(dims){};

	void set_dims(std::vector< size_t > const dims){ this->dims_ = dims;};
	std::vector< size_t > get_dims( void ){ return this->dims_;};

	
	std::vector< std::vector< int > > get_all_combinations( void )
	{
		return this->all_combinations_;
	};
	
	void compute_all_combinations( void )
	{
		size_t linear_range = 1;
		for( int direction=0; direction<dims_.size(); direction++ )
			linear_range *= dims_[direction];

		for( size_t i=0; i<linear_range; i++)
		{
			std::vector< int > curr_combination;

			for( int j=0; j<this->dims_.size(); j++)
			{
				size_t curr_idx = i;

				for(int k=0;k<j;k++)
					curr_idx = this->recurse(curr_idx, k, curr_combination);

				curr_idx = curr_idx % this->dims_[j];
				curr_combination.push_back(curr_idx);
			}
			this->all_combinations_.push_back( curr_combination );
		}
	};

private:

	std::vector< size_t > dims_;
	std::vector< std::vector< int > > all_combinations_;

	int recurse(size_t const idx, int const num_iter, std::vector< int > const curr_combination)
	{	
		int l = (idx - curr_combination[num_iter])/this->dims_[num_iter];
		return l;
	};

};