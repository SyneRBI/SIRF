/* file describing tissue parameter structs. 
these data containers should be filled by an xml parser.

author		johannes mayer
date		20180321
institute	PTB Berlin

*/


struct MRTissueParameter {

	float t1_miliseconds_;
	float t2_miliseconds_;
	float cs_ppm_;
};



struct PETTissueParameter {

	float attenuation_1_by_mm_;
	float suv_;
};


struct TissueParameter {

	int label_;
	char* name_;

	MRTissueParameter mr_tissue_;
	PETTissueParameter pet_tissue_;
};








