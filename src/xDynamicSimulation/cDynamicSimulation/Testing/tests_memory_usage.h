/* ================================================

Author: Johannes Mayer
Date: 2018.11.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once



#include <memory>
#include <vector>
#include <stdexcept>

#include "sirf/cGadgetron/gadgetron_data_containers.h"
#include <ismrmrd/ismrmrd.h>


void wait_for_time ( int const wait_time_s);


namespace tests_memory{


bool test_acquisition_memory( void );
bool test_acquisition_vector_memory( void );
bool test_downsizing_acquisition_memory( void );
bool test_acquisition_vector_ordering_memory( void );

bool test_ndarray_memory_managment( void );
bool tests_resizing_acquisition_memory( void );

bool tests_VD_h5_file_content( void );
}