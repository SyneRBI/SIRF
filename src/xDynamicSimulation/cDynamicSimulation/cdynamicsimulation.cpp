/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)

Copyright 2021 Physikalisch-Technische Bundesanstalt (PTB)

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

#include "sirf/cDynamicSimulation/cdynamicsimulation.h"

#include <memory>

#include "sirf/iUtilities/DataHandle.h"
#include "sirf/cDynamicSimulation/dynamicsimulation_x.h"
#include "sirf/cDynamicSimulation/contrastgenerator.h"

// using namespace sirf;


// extern "C"
// void*
// cGT_ISMRMRDAcquisitionsFromFile(const char* file)
// {
// 	if (!file_exists(file))
// 		return fileNotFound(file, __FILE__, __LINE__);
// 	try {
// 		shared_ptr<MRAcquisitionData>
// 			acquisitions(new AcquisitionsVector);
// 		acquisitions->read(file);
// 		return newObjectHandle<MRAcquisitionData>(acquisitions);
// 	}
// 	CATCH;
// }

extern "C"
void* cDS_MRDynamicSimulation(void* ptr_labels, const char* fname_xml)
{
	try {

        CAST_PTR(DataHandle, h_labels, ptr_labels);
		std::shared_ptr<LabelVolume> sptr_labels;
		getObjectSptrFromHandle<LabelVolume>(h_labels, sptr_labels);
		
        MRContrastGenerator cont_gen(*sptr_labels, fname_xml);
        
        std::shared_ptr<MRDynamicSimulation> 
			sptr_ds(new MRDynamicSimulation(cont_gen));

		return newObjectHandle<MRDynamicSimulation>(sptr_ds);

	}
	CATCH;
}