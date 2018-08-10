/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

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

#include <exception>
#include <iostream>
#include <string>

#include <boost/filesystem/operations.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/asio.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

using boost::asio::ip::tcp;

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>

#include "cgadgetron_shared_ptr.h"
#include "data_handle.h"
#include "gadgetron_data_containers.h"
#include "gadgetron_client.h"
//#include "iutilities.h" // causes problems with Matlab (cf. the same message below)
#include "cgadgetron_p.h"
#include "gadgetron_x.h"
#include "gadget_lib.h"
#include "chain_lib.h"

using namespace gadgetron;
using namespace sirf;

#define GRAB 1

#define NEW_GADGET(G) if (boost::iequals(name, G::class_name())) \
return newObjectHandle<G>();
#define NEW_GADGET_CHAIN(C) if (boost::iequals(name, C::class_name())) \
return newObjectHandle<C>();

shared_ptr<boost::mutex> Mutex::sptr_mutex_;

static void*
unknownObject(const char* obj, const char* name, const char* file, int line)
{
	DataHandle* handle = new DataHandle;
	std::string error = "Unknown ";
	error += obj;
	error += " '";
	error += name;
	error += "'";
	ExecutionStatus status(error.c_str(), file, line);
	handle->set(0, &status);
	return (void*)handle;
}

static void*
parameterNotFound(const char* name, const char* file, int line)
{
	DataHandle* handle = new DataHandle;
	std::string error = "Parameter ";
	error += name;
	error += " not found";
	ExecutionStatus status(error.c_str(), file, line);
	handle->set(0, &status);
	return (void*)handle;
}

static void*
fileNotFound(const char* name, const char* file, int line)
{
	DataHandle* handle = new DataHandle;
	std::string error = "File ";
	error += name;
	error += " not found";
	ExecutionStatus status(error.c_str(), file, line);
	handle->set(0, &status);
	return (void*)handle;
}

extern "C"
void* cGT_newObject(const char* name)
{
	try {
		if (boost::iequals(name, "Mutex"))
			return newObjectHandle<Mutex>();
		if (boost::iequals(name, "GTConnector"))
			return newObjectHandle<GTConnector>();
		if (boost::iequals(name, "CoilImages"))
			return newObjectHandle<CoilImagesVector>();
		if (boost::iequals(name, "AcquisitionModel"))
			return newObjectHandle<MRAcquisitionModel>();
		NEW_GADGET_CHAIN(GadgetChain);
		NEW_GADGET_CHAIN(AcquisitionsProcessor);
		NEW_GADGET_CHAIN(ImagesReconstructor);
		NEW_GADGET_CHAIN(ImagesProcessor);
		NEW_GADGET_CHAIN(RemoveOversamplingProcessor);
		NEW_GADGET_CHAIN(ExtractRealImagesProcessor);
		NEW_GADGET_CHAIN(SimpleReconstructionProcessor);
		NEW_GADGET_CHAIN(SimpleGRAPPAReconstructionProcessor);
		NEW_GADGET(IsmrmrdAcqMsgReader);
		NEW_GADGET(IsmrmrdAcqMsgWriter);
		NEW_GADGET(IsmrmrdImgMsgReader);
		NEW_GADGET(IsmrmrdImgMsgWriter);
		NEW_GADGET(NoiseAdjustGadget);
		NEW_GADGET(AsymmetricEchoAdjustROGadget);
		NEW_GADGET(RemoveROOversamplingGadget);
		NEW_GADGET(AcquisitionAccumulateTriggerGadget);
		NEW_GADGET(BucketToBufferGadget);
		NEW_GADGET(GenericReconCartesianReferencePrepGadget);
		NEW_GADGET(GenericReconCartesianGrappaGadget);
		NEW_GADGET(SimpleReconGadget);
		NEW_GADGET(GenericReconFieldOfViewAdjustmentGadget);
		NEW_GADGET(GenericReconImageArrayScalingGadget);
		NEW_GADGET(ImageArraySplitGadget);
		NEW_GADGET(ExtractGadget);
		NEW_GADGET(ComplexToFloatGadget);
		NEW_GADGET(FloatToShortGadget);
		NEW_GADGET(ImageFinishGadget);
		NEW_GADGET(AcquisitionFinishGadget);
		NEW_GADGET(SimpleReconGadgetSet);
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void* 
cGT_parameter(void* ptr, const char* obj, const char* name)
{
	try {
		if (boost::iequals(obj, "image"))
			return cGT_imageParameter(ptr, name);
		if (boost::iequals(obj, "acquisition"))
			return cGT_acquisitionParameter(ptr, name);
		if (boost::iequals(obj, "acquisitions"))
			return cGT_acquisitionsParameter(ptr, name);
		if (boost::iequals(obj, "gadget_chain")) {
			GadgetChain& gc = objectFromHandle<GadgetChain>(ptr);
			shared_ptr<aGadget> sptr = gc.gadget_sptr(name);
			if (sptr.get())
				return newObjectHandle(sptr);
			else {
				DataHandle* handle = new DataHandle;
				std::string error = "Gadget ";
				error += name;
				error += " not in the chain";
				ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
				handle->set(0, &status);
				return (void*)handle;
			}
		}
		if (boost::iequals(obj, "gadget")) {
			aGadget& g = objectFromHandle<aGadget>(ptr);
			std::string value = g.value_of(name);
			return charDataHandleFromCharData(value.c_str());
		}
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void*
cGT_setParameter(void* ptr, const char* obj, const char* par, const void* val)
{
	try {
		if (boost::iequals(obj, "coil_sensitivity"))
			return cGT_setCSParameter(ptr, par, val);
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void*
cGT_CoilSensitivities(const char* file)
{
	try {
		if (std::strlen(file) > 0) {
			shared_ptr<CoilSensitivitiesContainer>
				csms(new CoilSensitivitiesAsImages(file));
			return newObjectHandle<CoilSensitivitiesContainer>(csms);
		}
		else {
			shared_ptr<CoilSensitivitiesContainer>
				csms(new CoilSensitivitiesAsImages());
			return newObjectHandle<CoilSensitivitiesContainer>(csms);
		}
	}
	CATCH;
}

extern "C"
void*
cGT_setCSParameter(void* ptr, const char* par, const void* val)
{
	CAST_PTR(DataHandle, h_csms, ptr);
	CoilSensitivitiesContainer& csms =
		objectFromHandle<CoilSensitivitiesContainer>(h_csms);
	if (boost::iequals(par, "smoothness"))
		csms.set_csm_smoothness(dataFromHandle<int>(val));
	//csms.set_csm_smoothness(intDataFromHandle(val)); // causes problems with Matlab
	else
		return unknownObject("parameter", par, __FILE__, __LINE__);
	return new DataHandle;
}

extern "C"
void*
cGT_computeCoilImages(void* ptr_cis, void* ptr_acqs)
{
	try {
		CAST_PTR(DataHandle, h_csms, ptr_cis);
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		CoilImagesContainer& cis =
			objectFromHandle<CoilImagesContainer>(h_csms);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		cis.compute(acqs);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_computeCSMsFromCIs(void* ptr_csms, void* ptr_cis)
{
	try {
		CAST_PTR(DataHandle, h_csms, ptr_csms);
		CAST_PTR(DataHandle, h_cis, ptr_cis);
		CoilSensitivitiesContainer& csms =
			objectFromHandle<CoilSensitivitiesContainer>(h_csms);
		CoilImagesContainer& cis =
			objectFromHandle<CoilImagesContainer>(h_cis);
		csms.compute(cis);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_computeCoilSensitivities(void* ptr_csms, void* ptr_acqs)
{
	try {
		CAST_PTR(DataHandle, h_csms, ptr_csms);
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		CoilSensitivitiesContainer& csms =
			objectFromHandle<CoilSensitivitiesContainer>(h_csms);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		csms.compute(acqs);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_appendCSM
(void* ptr_csms, int nx, int ny, int nz, int nc, size_t ptr_re, size_t ptr_im)
{
	try {
		CAST_PTR(DataHandle, h_csms, ptr_csms);
		float* re = (float*)ptr_re;
		float* im = (float*)ptr_im;
		CoilSensitivitiesContainer& list =
			objectFromHandle<CoilSensitivitiesContainer>(h_csms);
		list.append_csm(nx, ny, nz, nc, re, im);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void
cGT_getCoilDataDimensions(void* ptr_csms, int csm_num, size_t ptr_dim)
{
	int* dim = (int*)ptr_dim;
	CAST_PTR(DataHandle, h_csms, ptr_csms);
	CoilDataContainer& list =
		objectFromHandle<CoilDataContainer>(h_csms);
	list.get_dim(csm_num, dim);
}

extern "C"
void
cGT_getCoilData(void* ptr_csms, int csm_num, size_t ptr_re, size_t ptr_im)
{
	float* re = (float*)ptr_re;
	float* im = (float*)ptr_im;
	CAST_PTR(DataHandle, h_csms, ptr_csms);
	CoilDataContainer& list =
		objectFromHandle<CoilDataContainer>(h_csms);
	list.get_data(csm_num, re, im);
}

extern "C"
void
cGT_getCoilDataAbs(void* ptr_csms, int csm_num, size_t ptr)
{
	float* v = (float*)ptr;
	CAST_PTR(DataHandle, h_csms, ptr_csms);
	CoilDataContainer& list =
		objectFromHandle<CoilDataContainer>(h_csms);
	list.get_data_abs(csm_num, v);
}

extern "C"
void*
cGT_AcquisitionModel(const void* ptr_acqs, const void* ptr_imgs)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		CAST_PTR(DataHandle, h_imgs, ptr_imgs);
		shared_ptr<MRAcquisitionData> acqs =
			objectSptrFromHandle<MRAcquisitionData>(h_acqs);
		shared_ptr<MRImageData> imgs =
			objectSptrFromHandle<MRImageData>(h_imgs);
		shared_ptr<MRAcquisitionModel> am(new MRAcquisitionModel(acqs, imgs));
		return newObjectHandle<MRAcquisitionModel>(am);
	}
	CATCH;
}

extern "C"
void*
cGT_setUpAcquisitionModel
(void* ptr_am, const void* ptr_acqs, const void* ptr_imgs)
{
	try {
		CAST_PTR(DataHandle, h_am, ptr_am);
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		CAST_PTR(DataHandle, h_imgs, ptr_imgs);
		MRAcquisitionModel& am = objectFromHandle<MRAcquisitionModel>(h_am);
		shared_ptr<MRAcquisitionData> sptr_acqs =
			objectSptrFromHandle<MRAcquisitionData>(h_acqs);
		shared_ptr<MRImageData> sptr_imgs =
			objectSptrFromHandle<MRImageData>(h_imgs);
		am.set_up(sptr_acqs, sptr_imgs);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_setAcquisitionModelParameter
(void* ptr_am, const char* name, const void* ptr)
{
	try {
		CAST_PTR(DataHandle, h_am, ptr_am);
		if (boost::iequals(name, "acquisition_template")) {
			CAST_PTR(DataHandle, handle, ptr);
			MRAcquisitionModel& am = objectFromHandle<MRAcquisitionModel>(h_am);
			shared_ptr<MRAcquisitionData> sptr =
				objectSptrFromHandle<MRAcquisitionData>(handle);
			am.set_acquisition_template(sptr);
		}
		else if (boost::iequals(name, "image_template")) {
			CAST_PTR(DataHandle, handle, ptr);
			MRAcquisitionModel& am = objectFromHandle<MRAcquisitionModel>(h_am);
			shared_ptr<MRImageData> sptr =
				objectSptrFromHandle<MRImageData>(handle);
			am.set_image_template(sptr);
		}
		else if (boost::iequals(name, "coil_sensitivity_maps")) {
			CAST_PTR(DataHandle, handle, ptr);
			MRAcquisitionModel& am = objectFromHandle<MRAcquisitionModel>(h_am);
			shared_ptr<CoilSensitivitiesContainer> sptr =
				objectSptrFromHandle<CoilSensitivitiesContainer>(handle);
			am.setCSMs(sptr);
		}
		else
			return unknownObject("parameter", name, __FILE__, __LINE__);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_setCSMs(void* ptr_am, const void* ptr_csms)
{
	try {
		CAST_PTR(DataHandle, h_am, ptr_am);
		CAST_PTR(DataHandle, h_csms, ptr_csms);
		MRAcquisitionModel& am = objectFromHandle<MRAcquisitionModel>(h_am);
		shared_ptr<CoilSensitivitiesContainer> sptr_csms =
			objectSptrFromHandle<CoilSensitivitiesContainer>(h_csms);
		am.setCSMs(sptr_csms);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs)
{
	try {
		CAST_PTR(DataHandle, h_am, ptr_am);
		CAST_PTR(DataHandle, h_imgs, ptr_imgs);
		MRAcquisitionModel& am = objectFromHandle<MRAcquisitionModel>(h_am);
		MRImageData& imgs = objectFromHandle<MRImageData>(h_imgs);
		shared_ptr<MRAcquisitionData> sptr_acqs = am.fwd(imgs);
		return newObjectHandle<MRAcquisitionData>(sptr_acqs);
	}
	CATCH;
}

extern "C"
void*
cGT_AcquisitionModelBackward(void* ptr_am, const void* ptr_acqs)
{
	try {
		CAST_PTR(DataHandle, h_am, ptr_am);
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		MRAcquisitionModel& am = objectFromHandle<MRAcquisitionModel>(h_am);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		shared_ptr<MRImageData> sptr_imgs = am.bwd(acqs);
		return newObjectHandle<MRImageData>(sptr_imgs);
	}
	CATCH;
}

extern "C"
void*
cGT_setAcquisitionsStorageScheme(const char* scheme)
{
	try{
		if (scheme[0] == 'f' || strcmp(scheme, "default") == 0)
			AcquisitionsFile::set_as_template();
		else
			AcquisitionsVector::set_as_template();
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_getAcquisitionsStorageScheme()
{
	return charDataHandleFromCharData
		(MRAcquisitionData::storage_scheme().c_str());
}

extern "C"
void*
cGT_orderAcquisitions(void* ptr_acqs)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		acqs.order();
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_ISMRMRDAcquisitionsFromFile(const char* file)
{
	if (!boost::filesystem::exists(file))
		return fileNotFound(file, __FILE__, __LINE__);
	try {
		shared_ptr<MRAcquisitionData> 
			acquisitions(new AcquisitionsFile(file));
		return newObjectHandle<MRAcquisitionData>(acquisitions);
	}
	CATCH;
}

extern "C"
void*
cGT_ISMRMRDAcquisitionsFile(const char* file)
{
	try {
		shared_ptr<MRAcquisitionData> 
			acquisitions(new AcquisitionsFile(file, true));
		return newObjectHandle<MRAcquisitionData>(acquisitions);
	}
	CATCH;
}

extern "C"
void*
cGT_processAcquisitions(void* ptr_proc, void* ptr_input)
{
	try {
		CAST_PTR(DataHandle, h_proc, ptr_proc);
		CAST_PTR(DataHandle, h_input, ptr_input);
		AcquisitionsProcessor& proc =
			objectFromHandle<AcquisitionsProcessor>(h_proc);
		MRAcquisitionData& input =
			objectFromHandle<MRAcquisitionData>(h_input);
		proc.process(input);
		shared_ptr<MRAcquisitionData> sptr_ac = proc.get_output();
		return newObjectHandle<MRAcquisitionData>(sptr_ac);
	}
	CATCH;
}

extern "C"
void*
cGT_acquisitionFromContainer(void* ptr_acqs, unsigned int acq_num)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		shared_ptr<ISMRMRD::Acquisition>
			sptr_acq(new ISMRMRD::Acquisition);
		acqs.get_acquisition(acq_num, *sptr_acq);
		return newObjectHandle<ISMRMRD::Acquisition>(sptr_acq);
	}
	CATCH;
}

extern "C"
void*
cGT_getAcquisitionsDimensions(void* ptr_acqs, size_t ptr_dim)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		shared_ptr<ISMRMRD::Acquisition>
			sptr_acq(new ISMRMRD::Acquisition);
		int num_reg_dim = acqs.get_acquisitions_dimensions(ptr_dim);
		return dataHandle(num_reg_dim);
	}
	CATCH;
}

extern "C"
void*
cGT_getAcquisitionsFlags(void* ptr_acqs, unsigned int n, size_t ptr_f)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		shared_ptr<ISMRMRD::Acquisition>
			sptr_acq(new ISMRMRD::Acquisition);
		acqs.get_acquisitions_flags(n, (int*)ptr_f);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_getAcquisitionsData
(void* ptr_acqs, unsigned int slice, size_t ptr_re, size_t ptr_im)
{
	try {
		float* re = (float*)ptr_re;
		float* im = (float*)ptr_im;
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		int n = acqs.get_acquisitions_data(slice, re, im);
		return dataHandle(n);
	}
	CATCH;
}

extern "C"
void*
cGT_setAcquisitionsData
(void* ptr_acqs, unsigned int na, unsigned int nc, unsigned int ns,
size_t ptr_re, size_t ptr_im)
{
	try {
		float* re = (float*)ptr_re;
		float* im = (float*)ptr_im;
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		int err = acqs.set_acquisition_data(na, nc, ns, re, im);
		DataHandle* handle = new DataHandle;
		if (err)
			handle->set_status("Mismatching acquisition dimensions", 
				__FILE__, __LINE__);
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void*
cGT_writeAcquisitions(void* ptr_acqs, const char* filename)
{
	try {
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(ptr_acqs);
		acqs.write(filename);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_acquisitionParameter(void* ptr_acq, const char* name)
{
	CAST_PTR(DataHandle, h_acq, ptr_acq);
	ISMRMRD::Acquisition& acq =
		objectFromHandle<ISMRMRD::Acquisition>(h_acq);
	if (boost::iequals(name, "version"))
		return dataHandle((int)acq.version());
	if (boost::iequals(name, "flags"))
		return dataHandle((int)acq.flags());
	if (boost::iequals(name, "measurement_uid"))
		return dataHandle((int)acq.measurement_uid());
	if (boost::iequals(name, "scan_counter"))
		return dataHandle((int)acq.scan_counter());
	if (boost::iequals(name, "acquisition_time_stamp"))
		return dataHandle((int)acq.acquisition_time_stamp());
	if (boost::iequals(name, "number_of_samples"))
		return dataHandle((int)acq.number_of_samples());
	if (boost::iequals(name, "available_channels"))
		return dataHandle((int)acq.available_channels());
	if (boost::iequals(name, "active_channels"))
		return dataHandle((int)acq.active_channels());
	if (boost::iequals(name, "discard_pre"))
		return dataHandle((int)acq.discard_pre());
	if (boost::iequals(name, "discard_post"))
		return dataHandle((int)acq.discard_post());
	if (boost::iequals(name, "center_sample"))
		return dataHandle((int)acq.center_sample());
	if (boost::iequals(name, "encoding_space_ref"))
		return dataHandle((int)acq.encoding_space_ref());
	if (boost::iequals(name, "trajectory_dimensions"))
		return dataHandle((int)acq.trajectory_dimensions());
	if (boost::iequals(name, "idx_kspace_encode_step_1"))
		return dataHandle((int)acq.idx().kspace_encode_step_1);
	if (boost::iequals(name, "idx_kspace_encode_step_2"))
		return dataHandle((int)acq.idx().kspace_encode_step_2);
	if (boost::iequals(name, "idx_average"))
		return dataHandle((int)acq.idx().average);
	if (boost::iequals(name, "idx_slice"))
		return dataHandle((int)acq.idx().slice);
	if (boost::iequals(name, "idx_contrast"))
		return dataHandle((int)acq.idx().contrast);
	if (boost::iequals(name, "idx_phase"))
		return dataHandle((int)acq.idx().phase);
	if (boost::iequals(name, "idx_repetition"))
		return dataHandle((int)acq.idx().repetition);
	if (boost::iequals(name, "idx_set"))
		return dataHandle((int)acq.idx().set);
	if (boost::iequals(name, "idx_segment"))
		return dataHandle((int)acq.idx().segment);
	if (boost::iequals(name, "physiology_time_stamp"))
		return dataHandle(acq.physiology_time_stamp());
	if (boost::iequals(name, "channel_mask"))
		return dataHandle(acq.channel_mask());
	if (boost::iequals(name, "sample_time_us"))
		return dataHandle((float)acq.sample_time_us());
	if (boost::iequals(name, "position"))
		return dataHandle((float*)acq.position());
	if (boost::iequals(name, "read_dir"))
		return dataHandle((float*)acq.read_dir());
	if (boost::iequals(name, "phase_dir"))
		return dataHandle((float*)acq.phase_dir());
	if (boost::iequals(name, "slice_dir"))
		return dataHandle((float*)acq.slice_dir());
	if (boost::iequals(name, "patient_table_position"))
		return dataHandle((float*)acq.patient_table_position());
	return parameterNotFound(name, __FILE__, __LINE__);
}

extern "C"
void*
cGT_acquisitionsParameter(void* ptr_acqs, const char* name)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		MRAcquisitionData& acqs =
			objectFromHandle<MRAcquisitionData>(h_acqs);
		if (boost::iequals(name, "undersampled"))
			return dataHandle((int)acqs.undersampled());
		return parameterNotFound(name, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void*
cGT_imageParameter(void* ptr_im, const char* name)
{
	try {
		ImageWrap& im = objectFromHandle<ImageWrap>(ptr_im);
		ISMRMRD::ImageHeader& head = im.head();
		if (boost::iequals(name, "version"))
			return dataHandle((int)head.version);
		if (boost::iequals(name, "flags"))
			return dataHandle((int)head.flags);
		if (boost::iequals(name, "data_type"))
			return dataHandle((int)head.data_type);
		if (boost::iequals(name, "measurement_uid"))
			return dataHandle((int)head.measurement_uid);
		if (boost::iequals(name, "channels"))
			return dataHandle((int)head.channels);
		if (boost::iequals(name, "average"))
			return dataHandle((int)head.average);
		if (boost::iequals(name, "slice"))
			return dataHandle((int)head.slice);
		if (boost::iequals(name, "contrast"))
			return dataHandle((int)head.contrast);
		if (boost::iequals(name, "phase"))
			return dataHandle((int)head.phase);
		if (boost::iequals(name, "repetition"))
			return dataHandle((int)head.repetition);
		if (boost::iequals(name, "set"))
			return dataHandle((int)head.set);
		if (boost::iequals(name, "acquisition_time_stamp"))
			return dataHandle((int)head.acquisition_time_stamp);
		if (boost::iequals(name, "image_type"))
			return dataHandle((int)head.image_type);
		if (boost::iequals(name, "image_index"))
			return dataHandle((int)head.image_index);
		if (boost::iequals(name, "image_series_index"))
			return dataHandle((int)head.image_series_index);
		if (boost::iequals(name, "attribute_string_len"))
			return dataHandle((int)head.attribute_string_len);
		if (boost::iequals(name, "matrix_size"))
			return dataHandle(head.matrix_size);
		if (boost::iequals(name, "physiology_time_stamp"))
			return dataHandle(head.physiology_time_stamp);
		if (boost::iequals(name, "field_of_view"))
			return dataHandle((float*)head.field_of_view);
		if (boost::iequals(name, "position"))
			return dataHandle((float*)head.position);
		if (boost::iequals(name, "read_dir"))
			return dataHandle((float*)head.read_dir);
		if (boost::iequals(name, "phase_dir"))
			return dataHandle((float*)head.phase_dir);
		if (boost::iequals(name, "slice_dir"))
			return dataHandle((float*)head.slice_dir);
		if (boost::iequals(name, "patient_table_position"))
			return dataHandle((float*)head.patient_table_position);
		return parameterNotFound(name, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void*
cGT_reconstructImages(void* ptr_recon, void* ptr_input)
{
	try {
		CAST_PTR(DataHandle, h_recon, ptr_recon);
		CAST_PTR(DataHandle, h_input, ptr_input);
		ImagesReconstructor& recon = objectFromHandle<ImagesReconstructor>(h_recon);
		MRAcquisitionData& input = objectFromHandle<MRAcquisitionData>(h_input);
		recon.process(input);
		shared_ptr<MRImageData> sptr_img = recon.get_output();
		return newObjectHandle<MRImageData>(sptr_img);
	}
	CATCH;

}

extern "C"
void*
cGT_reconstructedImages(void* ptr_recon)
{
	try {
		CAST_PTR(DataHandle, h_recon, ptr_recon);
		ImagesReconstructor& recon = objectFromHandle<ImagesReconstructor>(h_recon);
		shared_ptr<MRImageData> sptr_img = recon.get_output();
		return newObjectHandle<MRImageData>(sptr_img);
	}
	CATCH;

}

extern "C"
void*
cGT_readImages(const char* file)
{
	try {
		shared_ptr<MRImageData> sptr_img(new ImagesVector);
		sptr_img->read(file);
		return newObjectHandle<MRImageData>(sptr_img);
	}
	CATCH;
}

extern "C"
void*
cGT_processImages(void* ptr_proc, void* ptr_input)
{
	try {
		CAST_PTR(DataHandle, h_proc, ptr_proc);
		CAST_PTR(DataHandle, h_input, ptr_input);
		ImagesProcessor& proc = objectFromHandle<ImagesProcessor>(h_proc);
		MRImageData& input = objectFromHandle<MRImageData>(h_input);
		proc.process(input);
		shared_ptr<MRImageData> sptr_img = proc.get_output();
		return newObjectHandle<MRImageData>(sptr_img);
	}
	CATCH;

}

extern "C"
void*
cGT_selectImages(void* ptr_input, const char* attr, const char* target)
{
	try {
		CAST_PTR(DataHandle, h_input, ptr_input);
		MRImageData& input = objectFromHandle<MRImageData>(h_input);
		shared_ptr<MRImageData> sptr_img = input.clone(attr, target);
		return newObjectHandle<MRImageData>(sptr_img);
	}
	CATCH;
}

extern "C"
void*
cGT_writeImages(void* ptr_imgs, const char* out_file, const char* out_group)
{
	try {
		CAST_PTR(DataHandle, h_imgs, ptr_imgs);
		MRImageData& list = objectFromHandle<MRImageData>(h_imgs);
		list.write(out_file, out_group);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_imageWrapFromContainer(void* ptr_imgs, unsigned int img_num)
{
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	MRImageData& images = objectFromHandle<MRImageData>(h_imgs);
	return newObjectHandle<ImageWrap>(images.sptr_image_wrap(img_num));
}

extern "C"
void
cGT_getImageDim(void* ptr_img, size_t ptr_dim)
{
	int* dim = (int*)ptr_dim;
	ImageWrap& image = objectFromHandle<ImageWrap>(ptr_img);
	image.get_dim(dim);
}

extern "C"
void*
cGT_imageType(const void* ptr_img)
{
	try {
		ImageWrap& image = objectFromHandle<ImageWrap>(ptr_img);
		return dataHandle(image.type());
	}
	CATCH;
}

extern "C"
void
cGT_getImageDataAsFloatArray(void* ptr_img, size_t ptr_data)
{
	float* data = (float*)ptr_data;
	ImageWrap& image = objectFromHandle<ImageWrap>(ptr_img);
	image.get_data(data);
}

extern "C"
void
cGT_getImageDataAsComplexArray(void* ptr_img, size_t ptr_re, size_t ptr_im)
{
	float* re = (float*)ptr_re;
	float* im = (float*)ptr_im;
	ImageWrap& image = objectFromHandle<ImageWrap>(ptr_img);
	image.get_cmplx_data(re, im);
}

extern "C"
void
cGT_getImageDimensions(void* ptr_imgs, int img_num, size_t ptr_dim)
{
	int* dim = (int*)ptr_dim;
	MRImageData& images = objectFromHandle<MRImageData>(ptr_imgs);
	images.sptr_image_wrap(img_num)->get_dim(dim);
	//CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	//MRImageData& list = objectFromHandle<MRImageData>(h_imgs);
	//list.get_image_dimensions(img_num, dim);
}

extern "C"
void
cGT_getImagesDataAsFloatArray(void* ptr_imgs, size_t ptr_data)
{
	float* data = (float*)ptr_data;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	MRImageData& list = objectFromHandle<MRImageData>(h_imgs);
	list.get_images_data_as_float_array(data);
}

extern "C"
void
cGT_getImagesDataAsComplexArray(void* ptr_imgs, size_t ptr_re, size_t ptr_im)
{
	float* re = (float*)ptr_re;
	float* im = (float*)ptr_im;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	MRImageData& list = objectFromHandle<MRImageData>(h_imgs);
	list.get_images_data_as_complex_array(re, im);
}

extern "C"
void*
cGT_setComplexImagesData(void* ptr_imgs, size_t ptr_re, size_t ptr_im)
{
	try {
		float* re = (float*)ptr_re;
		float* im = (float*)ptr_im;
		CAST_PTR(DataHandle, h_imgs, ptr_imgs);
		MRImageData& list = objectFromHandle<MRImageData>(h_imgs);
		list.set_complex_images_data(re, im);
	}
	CATCH;
	return (void*)new DataHandle;
}

extern "C"
void*
cGT_imageDataType(const void* ptr_x, int im_num)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		MRImageData& x = objectFromHandle<MRImageData>(h_x);
		return dataHandle(x.image_data_type(im_num));
	}
	CATCH;
}

extern "C"
void*
cGT_dataItems(const void* ptr_x)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		aDataContainer<complex_float_t>& x = 
			objectFromHandle<aDataContainer<complex_float_t> >(h_x);
		return dataHandle(x.items());
	}
	CATCH;
}

extern "C"
void*
cGT_norm(const void* ptr_x)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		aDataContainer<complex_float_t>& x = 
			objectFromHandle<aDataContainer<complex_float_t> >(h_x);
		return dataHandle(x.norm());
	}
	CATCH;
}

extern "C"
void*
cGT_dot(const void* ptr_x, const void* ptr_y)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		CAST_PTR(DataHandle, h_y, ptr_y);
		aDataContainer<complex_float_t>& x = 
			objectFromHandle<aDataContainer<complex_float_t> >(h_x);
		aDataContainer<complex_float_t>& y = 
			objectFromHandle<aDataContainer<complex_float_t> >(h_y);
		return dataHandle(x.dot(y));
	}
	CATCH;
}

extern "C"
void*
cGT_axpby(
float ar, float ai, const void* ptr_x,
float br, float bi, const void* ptr_y
){
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		CAST_PTR(DataHandle, h_y, ptr_y);
		aDataContainer<complex_float_t>& x = 
			objectFromHandle<aDataContainer<complex_float_t> >(h_x);
		aDataContainer<complex_float_t>& y = 
			objectFromHandle<aDataContainer<complex_float_t> >(h_y);
		shared_ptr<aDataContainer<complex_float_t> > 
			sptr_z(x.new_data_container());
		complex_float_t a(ar, ai);
		complex_float_t b(br, bi);
		sptr_z->axpby(a, x, b, y);
		return newObjectHandle<aDataContainer<complex_float_t> >(sptr_z);
	}
	CATCH;
}

extern "C"
void*
cGT_multiply(const void* ptr_x, const void* ptr_y)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		CAST_PTR(DataHandle, h_y, ptr_y);
		aDataContainer<complex_float_t>& x =
			objectFromHandle<aDataContainer<complex_float_t> >(h_x);
		aDataContainer<complex_float_t>& y =
			objectFromHandle<aDataContainer<complex_float_t> >(h_y);
		shared_ptr<aDataContainer<complex_float_t> >
			sptr_z(x.new_data_container());
		sptr_z->multiply(x, y);
		return newObjectHandle<aDataContainer<complex_float_t> >(sptr_z);
	}
	CATCH;
}

extern "C"
void*
cGT_divide(const void* ptr_x, const void* ptr_y)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		CAST_PTR(DataHandle, h_y, ptr_y);
		aDataContainer<complex_float_t>& x =
			objectFromHandle<aDataContainer<complex_float_t> >(h_x);
		aDataContainer<complex_float_t>& y =
			objectFromHandle<aDataContainer<complex_float_t> >(h_y);
		shared_ptr<aDataContainer<complex_float_t> >
			sptr_z(x.new_data_container());
		sptr_z->divide(x, y);
		return newObjectHandle<aDataContainer<complex_float_t> >(sptr_z);
	}
	CATCH;
}

extern "C"
void*
cGT_setConnectionTimeout(void* ptr_con, unsigned int timeout_ms)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		con.set_timeout(timeout_ms);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_addReader(void* ptr_gc, const char* id, const void* ptr_r)
{
	try {
		CAST_PTR(DataHandle, h_gc, ptr_gc);
		CAST_PTR(DataHandle, h_r, ptr_r);
		GadgetChain& gc = objectFromHandle<GadgetChain>(h_gc);
		shared_ptr<aGadget>& g = objectSptrFromHandle<aGadget>(h_r);
		gc.add_reader(id, g);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_addWriter(void* ptr_gc, const char* id, const void* ptr_w)
{
	try {
		CAST_PTR(DataHandle, h_gc, ptr_gc);
		CAST_PTR(DataHandle, h_w, ptr_w);
		GadgetChain& gc = objectFromHandle<GadgetChain>(h_gc);
		shared_ptr<aGadget>& g = objectSptrFromHandle<aGadget>(h_w);
		gc.add_writer(id, g);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_addGadget(void* ptr_gc, const char* id, const void* ptr_g)
{
	try {
		CAST_PTR(DataHandle, h_gc, ptr_gc);
		CAST_PTR(DataHandle, h_g, ptr_g);
		GadgetChain& gc = objectFromHandle<GadgetChain>(h_gc);
		shared_ptr<aGadget>& g = objectSptrFromHandle<aGadget>(h_g);
		gc.add_gadget(id, g);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_setGadgetProperty(void* ptr_g, const char* prop, const char* value)
{
	try {
		CAST_PTR(DataHandle, h_g, ptr_g);
		aGadget& g = objectFromHandle<aGadget>(h_g);
		//std::cout << g.name() << std::endl;
		g.set_property(prop, value);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_setGadgetProperties(void* ptr_g, const char* props)
{
	try {
		CAST_PTR(DataHandle, h_g, ptr_g);
		aGadget& g = objectFromHandle<aGadget>(h_g);
		std::string in(props);
		std::string prop;
		std::string value;
		size_t n = in.length();
		size_t i, j;
		i = 0;
		for (;;) {
			j = in.find_first_not_of(" \t\n\v\f\r", i);
			if (j == std::string::npos)
				break;
			i = j;
			j = in.find_first_of("= \t\n\v\f\r", i);
			if (j == std::string::npos)
				j = n;
			prop = in.substr(i, j - i);
			//std::cout << prop << '\n';
			i = j;
			i = in.find_first_not_of("= \t\n\v\f\r", i);
			j = in.find_first_of(", \t\n\v\f\r", i);
			if (j == std::string::npos)
				j = n;
			value = in.substr(i, j - i);
			//std::cout << value << '\n';
			g.set_property(prop.c_str(), value.c_str());
			if (j < n && in[j] == ',')
				j++;
			i = j;
		}
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_configGadgetChain(void* ptr_con, void* ptr_gc)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		CAST_PTR(DataHandle, h_gc, ptr_gc);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		GadgetChain& gc = objectFromHandle<GadgetChain>(h_gc);
		std::string config = gc.xml();
		//std::cout << config << std::endl;
		con.send_gadgetron_configuration_script(config);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_registerImagesReceiver(void* ptr_con, void* ptr_img)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		CAST_PTR(DataHandle, h_img, ptr_img);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		shared_ptr<MRImageData> sptr_images =
			objectSptrFromHandle<MRImageData>(h_img);
		con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageCollector(sptr_images)));
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_connect(void* ptr_con, const char* host, const char* port)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		con.connect(host, port);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_sendConfigScript(void* ptr_con, const char* config)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		con.send_gadgetron_configuration_script(config);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_sendConfigFile(void* ptr_con, const char* file)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		con.send_gadgetron_configuration_file(file);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_sendParameters(void* ptr_con, const void* ptr_par)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		CAST_PTR(DataHandle, h_par, ptr_par);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		std::string& par = objectFromHandle<std::string>(h_par);
		//std::cout << par << std::endl;
		con.send_gadgetron_parameters(par);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_sendParametersString(void* ptr_con, const char* par)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		con.send_gadgetron_parameters(par);
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void* 
cGT_sendAcquisitions(void* ptr_con, void* ptr_dat)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		CAST_PTR(DataHandle, h_dat, ptr_dat);
	
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		Mutex mutex;
		boost::mutex& mtx = mutex();
		ISMRMRD::Dataset& ismrmrd_dataset = 
			objectFromHandle<ISMRMRD::Dataset>(h_dat);
	
		uint32_t acquisitions = 0;
		{
			mtx.lock();
			acquisitions = ismrmrd_dataset.getNumberOfAcquisitions();
			mtx.unlock();
		}

		//std::cout << acquisitions << " acquisitions" << std::endl;
	
		ISMRMRD::Acquisition acq_tmp;
		for (uint32_t i = 0; i < acquisitions; i++) {
			{
				boost::mutex::scoped_lock scoped_lock(mtx);
				ismrmrd_dataset.readAcquisition(i, acq_tmp);
			}
			con.send_ismrmrd_acquisition(acq_tmp);
		}
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_sendImages(void* ptr_con, void* ptr_img)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		CAST_PTR(DataHandle, h_img, ptr_img);

		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		MRImageData& images = objectFromHandle<MRImageData>(h_img);
		for (unsigned int i = 0; i < images.number(); i++) {
			ImageWrap& iw = images.image_wrap(i);
			con.send_wrapped_image(iw);
		}
	}
	CATCH;

	return (void*)new DataHandle;
}

extern "C"
void*
cGT_disconnect(void* ptr_con)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		con.send_gadgetron_close();
		con.wait();
	}
	CATCH;

	return (void*)new DataHandle;
}

