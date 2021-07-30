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
#include "sirf/Gadgetron/gadgetron_data_containers.h"

#include "sirf/cDynamicSimulation/dynamicsimulation_x.h"
#include "sirf/cDynamicSimulation/contrastgenerator.h"

using namespace sirf;


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
void* cDS_simulateData(void* ptr_sim)
{
	try {

		CAST_PTR(DataHandle, h_sim, ptr_sim);			
		
		aDynamicSimulation& sim = objectFromHandle<aDynamicSimulation>(h_sim);
		sim.simulate_data();
		return new DataHandle;

	}
	CATCH;
}

extern "C"
void* cDS_writeSimulationResults(const void* ptr_sim, const char* fname_with_ext)
{
	try {

		CAST_PTR(DataHandle, h_sim, ptr_sim);			
		aDynamicSimulation& sim = objectFromHandle<aDynamicSimulation>(h_sim);
		sim.write_simulation_results(fname_with_ext);

		return new DataHandle;

	}
	CATCH;
}

extern "C"
void* cDS_MRDynamicSimulation(const void* ptr_labels, const char* fname_xml)
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

extern "C"
void* cDS_setAcquisitionTemplateData(void* ptr_sim, const void* ptr_acqs)
{
	try {

		CAST_PTR(DataHandle, h_sim, ptr_sim);			
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		
		MRAcquisitionData& ad = objectFromHandle<MRAcquisitionData>(h_acqs);
		MRDynamicSimulation& sim = objectFromHandle<MRDynamicSimulation>(h_sim);
		sim.set_template_acquisition_data(ad);
		return new DataHandle;

	}
	CATCH;
}

extern "C"
void* cDS_setCoilmaps(void* ptr_sim, const void* ptr_csm)
{
	try {

		CAST_PTR(DataHandle, h_sim, ptr_sim);			
		MRDynamicSimulation& sim = objectFromHandle<MRDynamicSimulation>(h_sim);

		CAST_PTR(DataHandle, h_csm, ptr_csm);
		std::shared_ptr<CoilSensitivitiesVector> sptr_csm;
		getObjectSptrFromHandle(h_csm, sptr_csm);

		sim.set_coilmaps(sptr_csm); 

		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cDS_setSNR(void* ptr_sim, float const SNR)
{
	try {

		CAST_PTR(DataHandle, h_sim, ptr_sim);			
		MRDynamicSimulation& sim = objectFromHandle<MRDynamicSimulation>(h_sim);
		sim.set_SNR(SNR);

		return new DataHandle;

	}

	CATCH;
}

extern "C"
void* cDS_setNoiseLabel(void* ptr_sim, int const label)
{
	try {

		CAST_PTR(DataHandle, h_sim, ptr_sim);			
		MRDynamicSimulation& sim = objectFromHandle<MRDynamicSimulation>(h_sim);
		sim.set_noise_label(label);

		return new DataHandle;

	}

	CATCH;
}



// Dynamics

extern "C"
void* cDS_setDynamicSignal(void* ptr_dyn, size_t signal)
{
	try {

		// CAST_PTR(DataHandle, h_sim, ptr_sim);			
		// MRDynamicSimulation& sim = objectFromHandle<MRDynamicSimulation>(h_sim);
		// sim.set_noise_label(label);

		return new DataHandle;

	}

	CATCH;
}

// Motion Dynamics

extern "C"
void* cDS_addDisplacementField(void* ptr_dyn, const void* ptr_dvf)
{
	try {

		CAST_PTR(DataHandle, h_dyn, ptr_dyn);			
		MotionDynamic& dyn = objectFromHandle<MotionDynamic>(h_dyn);
		
		CAST_PTR(DataHandle, h_dvf, ptr_dvf);
		MotionFieldType& dvf = objectFromHandle<MotionFieldType>(h_dvf);
		
		dyn.add_displacement_field(dvf);

		return new DataHandle;

	}

	CATCH;
}