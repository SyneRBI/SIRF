/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


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

	MRTissueParameter mr_tiss;
	PETTissueParameter pet_tiss;

	tissue_pars.mr_tissue_ = mr_tiss;
	tissue_pars.pet_tissue_ = pet_tiss;

	return true;
}



bool test_get_MRTissueParameter_from_ptree()
{

	using boost::property_tree::ptree;
	ptree pt;

	float const input_t1 = 1000;
	float const input_t2 = 2000;
	float const input_cs = 4.3;


	pt.put("TissueParameterList.TissueParameter.MRTissueParameter.t1_miliseconds", input_t1);
	pt.put("TissueParameterList.TissueParameter.MRTissueParameter.t2_miliseconds", input_t2);
	pt.put("TissueParameterList.TissueParameter.MRTissueParameter.cs_ppm", input_cs);


	MRTissueParameter mr_tiss;

	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("TissueParameterList") )
	{
		if( v.first == "TissueParameter")		
		{
			mr_tiss = get_mrtissueparameter_from_ptree(v.second);
		}
	}


	bool parameter_set_correct = (mr_tiss.t1_miliseconds_ == input_t1);
	parameter_set_correct *= (mr_tiss.t2_miliseconds_ == input_t2);
	parameter_set_correct *= (mr_tiss.cs_ppm_ == input_cs);
	
	return (parameter_set_correct);

}



bool test_read_TissueParameter_label_from_xml( std::string const xml_filepath )
{

	TissueParameterList tissueList = read_TissueParameters_from_xml(xml_filepath);

	TissueParameter firstParam = tissueList[0];

	std::string const input_name = "Liver";
	int const input_label = 1;

	if ( input_name.compare(firstParam.name_)  || (firstParam.label_ != input_label) )
		return false;
	else
		return true;
}