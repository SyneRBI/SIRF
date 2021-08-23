/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include "sirf/cDynamicSimulation/dynamics.h"


namespace test_binprocessor{
    bool test_get_set_bins( void );
}

namespace test_surrogateprocessor{
    bool test_linear_interpolate_signal( void );

}

namespace test_contrastprocessor{

bool test_mr_contrast_motion_dyn_get_num_simul_states( void );

}

namespace test_motionprocessor{

bool test_motion_dynamic_counter( void );
bool test_motion_dynamic_temp_folder_setup( void );
bool test_motion_dynamic_set_motion_fields(void);
bool test_motion_dynamic_prep_motion_fields( void );
// bool test_mvf_vs_pet_img_quarternions( void );
bool test_motion_dynamic_temp_interpolate_dvfs( void );
bool test_motion_dynamic_save_gt_deformations( void );

bool test_nonisotropic_mvf_resampling( void );

}

namespace test_dynamic{

bool test_is_in_bin( void );
bool test_intersect_mr_acquisition_data( void );

bool test_bin_mr_acquisitions(void);
bool test_bin_pet_time_interval( void );

}
