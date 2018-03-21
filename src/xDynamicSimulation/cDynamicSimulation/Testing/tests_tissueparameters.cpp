/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#include "tissueparameters.h"
#include "tests_tissueparameters.h"


bool test_allocate_MRTissueParameter_successful(void)
{
	MRTissueParameter mr_tissue_pars;
	mr_tissue_pars.t1_miliseconds_ = 1000;
	mr_tissue_pars.t2_miliseconds_ = 200;
	mr_tissue_pars.cs_ppm_ = 4.23;

	return true;
}

bool test_allocate_PETTissueParameter_successful(void)
{
	PETTissueParameter pet_tissue_pars;
	pet_tissue_pars.attenuation_1_by_mm_ = 0.01;
	pet_tissue_pars.suv_ = 15;

	return true;
}

bool test_allocate_TissueParameter_successful(void)
{
	TissueParameter tissue_pars;
	tissue_pars.label_ = 1;
	tissue_pars.name_ = "Liver";

	return true;
}
