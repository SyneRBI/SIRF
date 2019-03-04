/* ================================================

Author: Johannes Mayer
Date: 2018.07.20
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */
#pragma once

#include <string>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "sirf/cGadgetron/gadgetron_data_containers.h"


#include "sirf/cGadgetron/gadgetron_x.h"
#include "sirf/cSTIR/stir_x.h"

#include "sirf/cGadgetron/encoding.h"

#include "sirf/cDynamicSimulation/tissuelabelmapper.h"
#include "sirf/cDynamicSimulation/contrastgenerator.h"
#include "sirf/cDynamicSimulation/dynamics.h"
#include "sirf/cDynamicSimulation/dynsim_noisegenerator.h"
#include "sirf/cDynamicSimulation/volume_orientator.h"




#define IMG_DATA_TYPE 7 // from ismrmrd enum ISMRMRD_CXFLOAT = 7





class aDynamicSimulation {

public:

	aDynamicSimulation(){};
	~aDynamicSimulation(){};

	std::string get_filename_rawdata( void )
	{ 
		return this-> filename_rawdata_; 
	}

	virtual void set_filename_rawdata( std::string const filename_template_rawdata ) 
	{ 
		this->filename_rawdata_ = filename_template_rawdata; 
	}

	virtual void simulate_dynamics( void ) = 0;
	virtual void write_simulation_results( const std::string& filename_output_with_extension ) = 0;
	virtual void save_ground_truth_displacements() = 0;


	virtual void acquire_raw_data( void ) = 0;



protected:

	std::string filename_rawdata_;

};



class MRDynamicSimulation : public aDynamicSimulation {

public:

	MRDynamicSimulation( MRContrastGenerator mr_cont_gen) : mr_cont_gen_( mr_cont_gen ) 
	{
		this->sptr_trajectory_ = std::shared_ptr< sirf::CartesianTrajectoryContainer >( new sirf::CartesianTrajectoryContainer() );
		this->vol_orientator_.set_readout_direction( sirf::ro_dir_z);
	};

	virtual void set_filename_rawdata( std::string const filename_template_rawdata );
	void write_simulation_results( const std::string& filename_output_with_extension );

	void save_ground_truth_displacements( void );

	void add_dynamic( std::shared_ptr<MRMotionDynamic> sptr_motion_dyn);
	void add_dynamic( std::shared_ptr<MRContrastDynamic> sptr_contrast_dyn); 

	ISMRMRD::IsmrmrdHeader get_ismrmrd_header( void ){ return this->hdr_;};
	
	void set_all_source_acquisitions(MRDataType& acquisitions );
	void set_SNR(float const SNR);
	void set_noise_label(size_t const label);

	void simulate_statics( void );
	void simulate_dynamics( void );


	void extract_hdr_information( void );
	void set_trajectory( std::shared_ptr<sirf::aTrajectoryContainer> sptr_trajectory);

	void set_coilmaps( ISMRMRD::Image< complex_float_t >& coilmaps );


	virtual void acquire_raw_data( void );
	

private:

	std::vector< std::shared_ptr<MRMotionDynamic> > motion_dynamics_;
	std::vector< std::shared_ptr<MRContrastDynamic> > contrast_dynamics_;

	GaussianNoiseGenerator noise_generator_;
	sirf::aVolumeOrientator vol_orientator_;

	ISMRMRD::IsmrmrdHeader hdr_;
	ISMRMRD::Image< complex_float_t > coilmaps_;

	MRDataType all_source_acquisitions_;
	MRDataType source_acquisitions_;
	MRDataType target_acquisitions_;
	
	MRContrastGenerator mr_cont_gen_;
	sirf::MRAcquisitionModel acq_model_;

	std::shared_ptr<sirf::aTrajectoryContainer> sptr_trajectory_;

	void shift_time_start_to_zero( void );
	void simulate_motion_dynamics( void );
	void simulate_contrast_dynamics( void );
	void simulate_simultaneous_motion_contrast_dynamics( void );
	void set_noise_scaling( std::shared_ptr<sirf::aTrajectoryContainer> sptr_traj );
};


class PETDynamicSimulation : public aDynamicSimulation{

public:
	PETDynamicSimulation( PETContrastGenerator pet_cont_gen ):aDynamicSimulation(), pet_cont_gen_(pet_cont_gen){};
		
	void simulate_statics( void );
	
	void simulate_dynamics( void );
	void simulate_dynamics( size_t const total_scan_time );

	void set_template_acquisition_data( void );
	void set_template_image_data( const std::string& filename_header_with_ext );

	void set_output_filename_prefix( const std::string& output_filename_prefix_);

	void add_dynamic( std::shared_ptr<PETMotionDynamic> sptr_motion_dyn);
	void add_dynamic( std::shared_ptr<PETContrastDynamic> sptr_contrast_dyn); 

	virtual void acquire_raw_data( void );
	
	void add_noise( void );
	void add_noise( float const scaling_factor );

	void write_simulation_results( const std::string& filename_output_with_extension );
	void save_ground_truth_displacements( void );

private:


	std::string output_filename_prefix_;

	void simulate_motion_dynamics(size_t const total_scan_time );	

	std::vector< std::shared_ptr<PETMotionDynamic> > motion_dynamics_;
	std::vector< std::shared_ptr<PETContrastDynamic> > contrast_dynamics_;

	std::shared_ptr<PoissonNoiseGenerator> sptr_noise_generator_;

	sirf::STIRImageData get_reduced_pet_img_in_template_format( const sirf::STIRImageData& full_size_img );

	PETContrastGenerator pet_cont_gen_;
	sirf::STIRImageData template_image_data_;

	sirf::PETAcquisitionModelUsingMatrix acq_model_;

	sirf::PETAcquisitionDataInFile source_acquisitions_;
	std::shared_ptr<sirf::PETAcquisitionData> sptr_target_acquisitions_;
	
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

	DimensionsType dims_;
	std::vector< DimensionsType > all_combinations_;

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

	
	int recurse(size_t const idx, int const num_iter, DimensionsType const curr_combination)
	{	
		int l = (idx - curr_combination[num_iter])/this->dims_[num_iter];
		return l;
	};
};
