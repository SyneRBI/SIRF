/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "encoding.h"

#include <gadgetron/hoNDArray.h>

ISMRMRD::NDArray<complex_float_t> aFullySampledFFT::get_k_data( void )
{
	return this->k_data_;	
}




FullySampledCartesianFFT::FullySampledCartesianFFT(ISMRMRD::IsmrmrdHeader hdr):
aFullySampledFFT( hdr )
{
}
