/* ================================================

Author: Johannes Mayer
Date: 2018.07.20
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once



namespace test_lin_combi_gen{

bool test_get_all_combinations( void );

}


namespace tests_mr_dynsim
{

bool test_acquisitionsvector_memory_management( void );

bool test_constructor( void );

void test_extract_hdr_information( void );

bool test_simulate_contrast_dynamics( void );
bool test_simulate_motion_dynamics( void );
bool test_simulate_simultaneous_motion_contrast_dynamics( void );

bool test_simulate_rpe_acquisition( void );

}

namespace test_pet_dynsim
{

bool test_constructor( void );
bool set_template_acquisition_data( void );
bool test_simulate_dynamics( void );


}