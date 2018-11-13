/* file describing tissue parameter structs. 
these data containers should be filled by an xml parser.

author		johannes mayer
date		20180321
institute	PTB Berlin

*/

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>


#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>



struct MRTissueParameter {

	float spin_density_percentH2O_;
	float t1_miliseconds_;
	float t2_miliseconds_;
	float cs_ppm_;
	~MRTissueParameter(){};

	MRTissueParameter operator+ (const MRTissueParameter& mr_to_add) const
	{
		MRTissueParameter mr_tiss;
		mr_tiss.spin_density_percentH2O_ = this->spin_density_percentH2O_ + mr_to_add.spin_density_percentH2O_;
		mr_tiss.t1_miliseconds_ = this->t1_miliseconds_ + mr_to_add.t1_miliseconds_;
		mr_tiss.t2_miliseconds_ = this->t2_miliseconds_ + mr_to_add.t2_miliseconds_;
		mr_tiss.cs_ppm_ = this->cs_ppm_ + mr_to_add.cs_ppm_;
		
		return mr_tiss;
	};

};



struct PETTissueParameter {

	float attenuation_1_by_cm_;
	float activity_kBq_ml_;
	~PETTissueParameter(){};

	PETTissueParameter operator+ (const PETTissueParameter& pet_to_add) const
	{
		PETTissueParameter pet_tiss;
		pet_tiss.attenuation_1_by_cm_ = this->attenuation_1_by_cm_ + pet_to_add.attenuation_1_by_cm_;
		pet_tiss.activity_kBq_ml_ = this->activity_kBq_ml_ + pet_to_add.activity_kBq_ml_;
		return pet_tiss;
	}

};




struct TissueParameter {

	int label_;
	std::string name_;

	MRTissueParameter mr_tissue_;
	PETTissueParameter pet_tissue_;

	~TissueParameter(){};

	TissueParameter operator+ (const TissueParameter& tiss_to_add) const
	{
		TissueParameter tiss;
		tiss.label_ = this->label_;
		tiss.name_ = this->name_;
		
		tiss.mr_tissue_ = this->mr_tissue_ + tiss_to_add.mr_tissue_;
		tiss.pet_tissue_ = this->pet_tissue_ + tiss_to_add.pet_tissue_;

		return tiss;
	};
};


MRTissueParameter operator* (float const x, const MRTissueParameter& a_mr);
PETTissueParameter operator* (float const x, const PETTissueParameter& a_pet);
TissueParameter operator* (float const x, const TissueParameter& a_tiss);



typedef std::vector< TissueParameter > TissueParameterList;


TissueParameterList read_TissueParameters_from_xml(const std::string& xml_filepath);

MRTissueParameter get_mrtissueparameter_from_ptree(boost::property_tree::ptree const pt);
PETTissueParameter get_pettissueparameter_from_ptree(boost::property_tree::ptree const pt);


bool check_label_uniqueness( const TissueParameterList& tiss_list);