/* ================================================

Author: Johannes Mayer
Date: 2018.11.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once



#include <memory>
#include <vector>
#include <stdexcept>

#include "gadgetron_data_containers.h"
#include <ismrmrd/ismrmrd.h>


namespace tests_memory{


bool test_acquisition_memory( void );
bool test_acquisition_vector_memory( void );

}