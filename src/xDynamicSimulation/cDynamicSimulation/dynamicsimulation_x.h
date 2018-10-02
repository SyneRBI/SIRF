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

#include "encoding.h"

#include "tissuelabelmapper.h"
#include "contrastgenerator.h"
#include "dynamics.h"
#include "dynsim_noisegenerator.h"
#include "volume_orientator.h"




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

	void add_dynamic( const MotionDynamic& motion_dyn) 
	{
		this->motion_dynamics_.push_back(motion_dyn);
	};


	void add_dynamic( const ContrastDynamic& cont_dyn) 
	{
		this->contrast_dynamics_.push_back(cont_dyn);
	};

	virtual void acquire_raw_data( void ) = 0;


protected:

	std::vector< MotionDynamic > motion_dynamics_;
	std::vector< ContrastDynamic > contrast_dynamics_;

	std::string filename_rawdata_;

};


typedef sirf::AcquisitionsVector MRDataContainerType;

class MRDynamicSimulation : public aDynamicSimulation {

public:

	MRDynamicSimulation( MRContrastGenerator mr_cont_gen) : mr_cont_gen_( mr_cont_gen ) 
	{
		this->sptr_trajectory_ = std::shared_ptr< sirf::CartesianTrajectoryContainer >( new sirf::CartesianTrajectoryContainer() );
		this->vol_orientator_.set_readout_direction( sirf::ro_dir_z);
	};
	void write_simulation_results( std::string const filename_output_with_extension );


	ISMRMRD::IsmrmrdHeader get_ismrmrd_header( void ){ return this->hdr_;};
	
	void set_all_source_acquisitions(MRDataContainerType acquisitions );
	void set_noise_width(float const sigma);
	void set_SNR(float const SNR);

	void simulate_statics( void );
	void simulate_dynamics( void );


	void extract_hdr_information( void );
	void set_trajectory( std::shared_ptr<sirf::aTrajectoryContainer> sptr_trajectory);

	void set_coilmaps( ISMRMRD::Image< complex_float_t > &coilmaps );


	virtual void acquire_raw_data( void );

private:

	GaussianNoiseGenerator noise_generator_;
	sirf::aVolumeOrientator vol_orientator_;

	ISMRMRD::IsmrmrdHeader hdr_;
	ISMRMRD::Image< complex_float_t > coilmaps_;

	MRDataContainerType all_source_acquisitions_;
	MRDataContainerType source_acquisitions_;
	MRDataContainerType target_acquisitions_;
	
	MRContrastGenerator mr_cont_gen_;
	sirf::MRAcquisitionModel acq_model_;

	std::shared_ptr<sirf::aTrajectoryContainer> sptr_trajectory_;

	void simulate_motion_dynamics( void );
	void simulate_contrast_dynamics( void );
	void simulate_simultaneous_motion_contrast_dynamics( void );

};


class PETDynamicSimulation : public aDynamicSimulation{

public:
	PETDynamicSimulation( PETContrastGenerator pet_cont_gen ):aDynamicSimulation(), pet_cont_gen_(pet_cont_gen){};
		
	void simulate_dynamics( void );

	void set_template_acquisition_data( void );

	virtual void acquire_raw_data( void );

	void write_simulation_results( std::string const filename_output_with_extension );

	

private:

	PETContrastGenerator pet_cont_gen_;
	sirf::PETAcquisitionModelUsingMatrix acq_model_;

	sirf::PETAcquisitionDataInFile source_acquisitions_;
	std::shared_ptr<sirf::PETAcquisitionData> target_acquisitions_;
	
};







typedef std::vector<int> DimensionsType;

class LinearCombiGenerator{

public:
	LinearCombiGenerator(DimensionsType const dims): dims_(dims){};

	void set_dims(DimensionsType const dims){ this->dims_ = dims;};
	DimensionsType get_dims( void ){ return this->dims_;};

	
	std::vector< DimensionsType > get_all_combinations( void )
	{
		this->compute_all_combinations();
		return this->all_combinations_;
	};
	
	size_t get_num_total_combinations()
	{
		size_t linear_range = 1;
		for( int direction=0; direction<dims_.size(); direction++ )
			linear_range *= dims_[direction];

		return linear_range;
	}

private:

	void compute_all_combinations( void )
	{
		size_t const linear_range = this->get_num_total_combinations();

		for( size_t i=0; i<linear_range; i++)
		{
			DimensionsType curr_combination;

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

	DimensionsType dims_;
	std::vector< DimensionsType > all_combinations_;

	int recurse(size_t const idx, int const num_iter, DimensionsType const curr_combination)
	{	
		int l = (idx - curr_combination[num_iter])/this->dims_[num_iter];
		return l;
	};
};
