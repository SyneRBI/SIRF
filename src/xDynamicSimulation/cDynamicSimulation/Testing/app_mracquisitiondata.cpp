/* ================================================

Author: Johannes Mayer
Date: 2019.10.07
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "app_mracquisitiondata.h" 

#include "sirf/cGadgetron/gadgetron_data_containers.h"
#include <ismrmrd/ismrmrd.h>


void apps_johannesmayer::omitt_first_acquisition(std::string const fname_ismrmrd)
{

	size_t const first_acq = 1; //start with this acquisition

	sirf::AcquisitionsVector acq_vec;
	acq_vec.read( fname_ismrmrd );

	sirf::AcquisitionsVector coherent_acq_vec;
	coherent_acq_vec.copy_acquisitions_info(acq_vec);

	for(size_t i_acq=first_acq; i_acq<acq_vec.number(); i_acq++)
	{	
		ISMRMRD::Acquisition temp_acq;
	    acq_vec.get_acquisition(i_acq, temp_acq);
		coherent_acq_vec.append_acquisition(temp_acq);
	}


	std::string fname_output = fname_ismrmrd;
	
	std::string const substring = ".h5";
	std::string const replacementstring = "_coherent_readout.h5";
	
	fname_output.replace(fname_output.find(substring), substring.length(), replacementstring);
	coherent_acq_vec.write(fname_output);
	
}