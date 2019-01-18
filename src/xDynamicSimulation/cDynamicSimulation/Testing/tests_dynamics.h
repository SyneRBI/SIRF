/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include "sirf/cDynamicSimulation/dynamics.h"




namespace test_dynamic{

bool test_is_in_bin( void );
bool test_intersect_mr_acquisition_data( void );


bool test_get_set_bins( void );
bool test_linear_interpolate_signal( void );
bool test_bin_mr_acquisitions(void);

bool test_motion_dynamic_counter( void );
bool test_motion_dynamic_temp_folder_setup( void );
bool test_motion_dynamic_set_motion_fields(void);
bool test_motion_dynamic_prep_motion_fields( void );
bool test_mvf_vs_pet_img_quarternions( void );
bool test_motion_dynamic_temp_interpolate_dvfs( void );

bool test_mr_contrast_motion_dyn_get_num_simul_states( void );



bool test_bin_pet_time_interval( void );

	 
}
