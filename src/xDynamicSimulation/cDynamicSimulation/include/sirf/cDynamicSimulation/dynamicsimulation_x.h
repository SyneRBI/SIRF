/* ================================================

Author: Johannes Mayer
Date: 2018.07.20
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */
#pragma once

#include <string>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"
#include "sirf/STIR/stir_x.h"

#include "sirf/cDynamicSimulation/tissuelabelmapper.h"
#include "sirf/cDynamicSimulation/contrastgenerator.h"
#include "sirf/cDynamicSimulation/dynamics.h"
#include "sirf/cDynamicSimulation/dynsim_noisegenerator.h"
#include "sirf/cDynamicSimulation/dynsim_deformer.h"


#define IMG_DATA_TYPE 7 // from ismrmrd enum ISMRMRD_CXFLOAT = 7



typedef std::vector<int> DimensionsType;

/*!
	\brief A class to generate combinations of integers in a vector. 
*/
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

class aDynamicSimulation {

public:

	aDynamicSimulation(){};
	~aDynamicSimulation(){};

	virtual void simulate_data( void ) = 0;
	virtual void write_simulation_results( const std::string& filename_output_with_extension ) = 0;
	
	virtual void save_ground_truth_displacements(void) const = 0;
	virtual void acquire_raw_data( void ) = 0;
	
protected:
	DynamicSimulationDeformer dsd_;

};


class MRDynamicSimulation : public aDynamicSimulation {

public:

	MRDynamicSimulation(MRContrastGenerator mr_cont_gen) : mr_cont_gen_( mr_cont_gen ) 
	{};

	void write_simulation_results( const std::string& filename_output_with_extension );

	
	void set_acquisition_template_rawdata(const sirf::MRAcquisitionData& acquisitions);
	void set_contrast_template_rawdata(const sirf::MRAcquisitionData& acquisitions);

	void set_SNR(float const SNR);
	void set_noise_label(int const label);

	void add_dynamic( std::shared_ptr<MRMotionDynamic> sptr_motion_dyn){
		this->motion_dynamics_.push_back(sptr_motion_dyn);
	}
	void add_dynamic( std::shared_ptr<MRContrastDynamic> sptr_contrast_dyn){
		this->contrast_dynamics_.push_back(sptr_contrast_dyn);
	} 

	void add_dynamic( std::shared_ptr<ExternalMRContrastDynamic> sptr_ext_contrast_dyn){
		this->external_contrast_.push_back(sptr_ext_contrast_dyn);
	}

	void simulate_statics( void );
	void simulate_data( void );

	void set_coilmaps(const std::shared_ptr<sirf::CoilSensitivitiesVector> sptr_csm)
	{
		this->acq_model_.set_csm(sptr_csm);
	};

	virtual void acquire_raw_data( void );
	virtual void save_ground_truth_displacements() const;
	
	virtual void save_groud_truth_parameter_maps( const std::string prefix_output );

	virtual void set_offset_transformation(const sirf::AffineTransformation<float>& trafo)
	{
		dsd_.set_offset_transformation(trafo);
	}

	TissueParameter get_petmr_tissue_parameter(LabelType label){
		return mr_cont_gen_.get_petmr_tissue_parameter(label);
	}

private:

	std::vector< std::shared_ptr<MRMotionDynamic> > motion_dynamics_;
	std::vector< std::shared_ptr<MRContrastDynamic> > contrast_dynamics_;
	std::vector< std::shared_ptr<ExternalMRContrastDynamic> > external_contrast_;
	
	GaussianNoiseGenerator noise_generator_;

	std::shared_ptr<sirf::MRAcquisitionData> sptr_source_acquisitions_;
	std::shared_ptr<sirf::MRAcquisitionData> sptr_template_data_;
	std::shared_ptr<sirf::MRAcquisitionData> sptr_simul_data_;

	MRContrastGenerator mr_cont_gen_;
	sirf::MRAcquisitionModel acq_model_;

	sirf::AcquisitionsVector 
		get_acquisitions_for_motionstate(DimensionsType current_combination) const;

	std::vector<sirf::NiftiImageData3DDeformation<float> >
		get_motionfields_for_motionstate(DimensionsType current_combination) const;

	LinearCombiGenerator prepare_motion_information(void);

	void update_tissue_parameters(TimeAxisType current_time_point);


	void shift_time_start_to_zero( void );
	void simulate_motion_dynamics( void );
	void simulate_contrast_dynamics( void );
	void simulate_simultaneous_motion_contrast_dynamics( void );
	void simulate_external_contrast_motion_dynamics( void );

	void set_noise_scaling();
};


class PETDynamicSimulation : public aDynamicSimulation{

public:
	PETDynamicSimulation( PETContrastGenerator pet_cont_gen ):aDynamicSimulation(), pet_cont_gen_(pet_cont_gen){};
		
	void simulate_statics( void );
	
	void simulate_data( void );
	void simulate_data( size_t const total_scan_time );

	std::string get_filename_rawdata( void )
	{ 
		return this-> filename_rawdata_; 
	}

	virtual void set_filename_rawdata( std::string const filename_template_rawdata ) 
	{ 
		this->filename_rawdata_ = filename_template_rawdata; 
	}

	void set_template_acquisition_data( void );
	void set_template_image_data( const std::string& filename_header_with_ext );

	void set_output_filename_prefix( const std::string& output_filename_prefix_);

	virtual void acquire_raw_data( void );
	
	void add_dynamic( std::shared_ptr<PETMotionDynamic> sptr_motion_dyn){
		this->motion_dynamics_.push_back(sptr_motion_dyn);
	}
	void add_dynamic( std::shared_ptr<PETContrastDynamic> sptr_contrast_dyn){
		this->contrast_dynamics_.push_back(sptr_contrast_dyn);
	} 

	void add_noise( void );
	void add_noise( float const scaling_factor );

	void write_simulation_results( const std::string& filename_output_with_extension );
	virtual void save_ground_truth_displacements() const;
	virtual void save_groud_truth_parameter_maps() const
	{
		throw std::runtime_error("Not implemented yet");
	}

private:

	std::string filename_rawdata_;
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





