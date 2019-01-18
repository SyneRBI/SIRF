/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */
// file containing auxiliary functions

#pragma once

#include <string>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "sirf/cGadgetron/gadgetron_data_containers.h"
#include "sirf/cGadgetron/gadgetron_image_wrap.h"

#include "auxiliary_input_output.h"

#include "phantom_input.h"
#include "tissueparameters.h"
#include "tissuelabelmapper.h"
#include "contrastgenerator.h"
#include "dynamics.h"
#include "sirf/cGadgetron/encoding.h"


#include "test_input_filenames.h"

// volume sizes

#define MOCK_FOV 256
#define MOCK_DATA_MATRIX_SIZE 64
#define MOCK_DATA_NUM_CHANNELS 4
#define MOCK_DATA_RO_OVERSAMPLING 1
#define MOCK_IMAGE_TYPE 5 // from ismrmrd enum ISMRMRD_IMTYPE_COMPLEX   = 5
#define MOCK_DATA_TYPE 7 // from ismrmrd enum ISMRMRD_CXFLOAT = 7
#define MOCK_FIELD_STRENGTH 1

// mock signal
#define MOCK_NUM_SIG_PTS 10


namespace aux_test
{

	TissueParameterList get_mock_tissue_param_list( void );
	LabelArray get_mock_label_array( void );

	TissueParameter get_mock_tissue_parameter( void );
	PETTissueParameter get_mock_PET_tissue_parameter( void );
	MRTissueParameter get_mock_MR_tissue_parameter( void );

	std::pair< TissueParameter, TissueParameter> get_mock_contrast_signal_extremes( void );

	MRContrastGenerator get_mock_mr_contrast_generator( void );
	PETContrastGenerator get_mock_pet_contrast_generator( void );

	ISMRMRD::IsmrmrdHeader get_mock_ismrmrd_header( void );
	std::string get_serialized_mock_ismrmrd_header( void );


	ISMRMRD::AcquisitionSystemInformation get_mock_acquisition_system_information( void );
	ISMRMRD::SequenceParameters get_mock_sequence_parameters( void );
	ISMRMRD::ExperimentalConditions get_mock_experimental_conditions( void );
	std::vector< ISMRMRD::Encoding > get_mock_encoding_vector( void );

	ISMRMRD::ExperimentalConditions get_mock_experimental_conditions( void );
	ISMRMRD::EncodingSpace get_mock_encoded_space( void );
	ISMRMRD::EncodingSpace get_mock_recon_space( void );
	ISMRMRD::EncodingLimits get_mock_encoding_limits( void );



	ISMRMRD::NDArray<complex_float_t> get_mock_ndarray_with_cube( void );
	ISMRMRD::Image< complex_float_t > get_mock_ismrmrd_image_with_cube( void );
	ISMRMRD::Image< float > get_mock_ismrmrd_image_with_gradients( void );

	ISMRMRD::NDArray<complex_float_t> get_mock_csm( void );
	ISMRMRD::Image<complex_float_t> get_mock_gaussian_csm( std::vector<size_t> vol_dims, int const num_coils );
	sirf::CoilDataAsCFImage get_mock_coildata_as_cfimage( void );

	ISMRMRD::AcquisitionHeader get_mock_acquisition_header( void );	
	sirf::AcquisitionsVector get_mock_acquisition_vector ( ISMRMRD::IsmrmrdHeader );	

	sirf::RPETrajectoryContainer get_mock_radial_trajectory(size_t const NRad, size_t const NAng);



	SignalContainer get_generic_respiratory_signal( sirf::AcquisitionsVector &acq_vec);
	SignalContainer get_generic_cardiac_signal( sirf::AcquisitionsVector &acq_vec);
	SignalContainer get_generic_contrast_inflow_signal( sirf::AcquisitionsVector &acq_vec);
	SignalContainer get_generic_contrast_in_and_outflow_signal( sirf::AcquisitionsVector &acq_vec);



	SignalContainer get_mock_motion_signal( void );

	SignalContainer get_mock_sinus_signal( sirf::AcquisitionsVector &acq_vec, TimeAxisType const period_duration_ms);
	SignalContainer get_mock_sawtooth_signal( sirf::AcquisitionsVector acq_vec, TimeAxisType const period_duration_ms);

	

	template <typename T> bool equal_array_content( ISMRMRD::NDArray<T> one_array, ISMRMRD::NDArray<T> other_array)
	{
		size_t const num_elements = one_array.getNumberOfElements();
		size_t const num_elements_other = other_array.getNumberOfElements();

		if( num_elements != num_elements_other)
		{
			return false;
		}
		else
		{
			bool content_is_equal = true;
			for(int i=0; i<num_elements; i++)
			{
				content_is_equal *= ( one_array(i) == other_array(i) );
			}
			return content_is_equal;
		}
	};


}// namespace aux_test
