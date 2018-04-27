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


};



struct PETTissueParameter {

	float attenuation_1_by_mm_;
	float suv_;
	~PETTissueParameter(){};

};


struct TissueParameter {

	int label_;
	std::string name_;

	MRTissueParameter mr_tissue_;
	PETTissueParameter pet_tissue_;

	~TissueParameter(){};
};



typedef std::vector< TissueParameter > TissueParameterList;


TissueParameterList read_TissueParameters_from_xml(std::string const xml_filepath);

MRTissueParameter get_mrtissueparameter_from_ptree(boost::property_tree::ptree const pt);
PETTissueParameter get_pettissueparameter_from_ptree(boost::property_tree::ptree const pt);


bool check_label_uniqueness( TissueParameterList const tiss_list);