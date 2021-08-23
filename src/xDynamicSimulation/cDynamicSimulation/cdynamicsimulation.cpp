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
#include "sirf/cDynamicSimulation/dynamics.h"

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

extern "C"
void* cDS_addMRMotionDynamic(void* ptr_sim, void* ptr_dyn)
{
	try {

		CAST_PTR(DataHandle, h_sim, ptr_sim);			
		MRDynamicSimulation& sim = objectFromHandle<MRDynamicSimulation>(h_sim);
		
		CAST_PTR(DataHandle, h_dyn, ptr_dyn);			
		std::shared_ptr<MRMotionDynamic> sptr_dyn;
		getObjectSptrFromHandle(h_dyn, sptr_dyn);
		
		sim.add_dynamic(sptr_dyn);

		return new DataHandle;

	}

	CATCH;
}

// Dynamics
extern "C"
void* cDS_setDynamicSignal(void* ptr_dyn, size_t ptr_time, size_t ptr_signal, int const num_points)
{
	try {

		float* time = (float*) ptr_time;
		float* signal = (float*) ptr_signal;

		typedef std::pair<TimeAxisType, SignalAxisType> SignalPoint;
		typedef std::vector< SignalPoint > SignalContainer;

        SignalContainer surrogate(num_points);
        for(size_t i=0; i<surrogate.size(); ++i)
		{
            surrogate.at(i).first = *(time+i);
			surrogate.at(i).second = *(signal+i);
		}

		CAST_PTR(DataHandle, h_dyn, ptr_dyn);			
		Dynamic& dyn = objectFromHandle<Dynamic>(h_dyn);

		dyn.set_dynamic_signal(surrogate);

		return new DataHandle;
	}

	CATCH;
}

extern "C"
void* cDS_setCyclicality(void* ptr_dyn, bool const cyc)
{
	try {

		CAST_PTR(DataHandle, h_dyn, ptr_dyn);			
		Dynamic& dyn = objectFromHandle<Dynamic>(h_dyn);
		dyn.set_cyclicality(cyc);

		return new DataHandle;
	}

	CATCH;
}

// Motion Dynamics

extern "C"
void* cDS_addMRDisplacementField(void* ptr_dyn, const void* ptr_dvf)
{
	try {

		CAST_PTR(DataHandle, h_dyn, ptr_dyn);			
		MRMotionDynamic& dyn = objectFromHandle<MRMotionDynamic>(h_dyn);
		
		CAST_PTR(DataHandle, h_dvf, ptr_dvf);
		MotionFieldType& dvf = objectFromHandle<MotionFieldType>(h_dvf);
		
		dyn.add_displacement_field(dvf);

		return new DataHandle;

	}

	CATCH;
}

// MR Dynamics

extern "C"
void* cDS_setMRAcquisitions(void* ptr_dyn, void* ptr_ad)
{
	try {

		CAST_PTR(DataHandle, h_dyn, ptr_dyn);			
		MRDynamic& dyn = objectFromHandle<MRDynamic>(h_dyn);
		
		CAST_PTR(DataHandle, h_ad, ptr_ad);
		MRAcquisitionData& ad = objectFromHandle<MRAcquisitionData>(h_ad);

		dyn.bin_mr_acquisitions(ad);
		
		return new DataHandle;
		
	}
	CATCH;
}

extern "C"
void* cDS_MRMotionDynamic( int const num_states )
{
	try {

		std::shared_ptr<MRDynamic> 
			sptr_dyn(new MRMotionDynamic(num_states));

		return newObjectHandle<MRDynamic>(sptr_dyn);
	}

	CATCH;
}