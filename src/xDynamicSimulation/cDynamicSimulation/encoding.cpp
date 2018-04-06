/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "encoding.h"

#include <gadgetron/hoNDArray.h>


FullySampledCartesianFFT::FullySampledCartesianFFT(ISMRMRD::IsmrmrdHeader hdr):
aFullySampledFFT( hdr )
{
	Gadgetron::hoNDArray< float > some_arr;

}