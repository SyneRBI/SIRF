#include <exception>
#include <iostream>
#include <string>

#include <boost/filesystem/operations.hpp>
#include <boost/shared_ptr.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>

#include "gadgetron_data_containers.h"
#include "gadgetron_client.h"
#include "data_handle.h"
//#include "iutilities.h" // causes problems with Matlab (cf. the same message below)
#include "cgadgetron_par.h"
#include "xgadgetron.h"
#include "gadget_lib.h"
#include "chain_lib.h"

#define GRAB 1

#define NEW_GADGET(G) if (boost::iequals(name, G::class_name())) \
return newObjectHandle<aGadget, G>();
#define NEW_GADGET_CHAIN(C) if (boost::iequals(name, C::class_name())) \
return newObjectHandle<GadgetChain, C>();

boost::shared_ptr<boost::mutex> Mutex::sptr_mutex_;

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
			return newObjectHandle<Mutex, Mutex>();
		if (boost::iequals(name, "GTConnector"))
			return newObjectHandle<GTConnector, GTConnector>();
		//if (boost::iequals(name, "string"))
		//	return newObjectHandle<std::string, std::string>();
		if (boost::iequals(name, "ImagesList"))
			return newObjectHandle<ImagesContainer, ImagesList>();
		if (boost::iequals(name, "CoilImagesList"))
			return newObjectHandle<CoilImagesContainer, CoilImagesList>();
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
		if (boost::iequals(obj, "acquisition"))
			return cGT_acquisitionParameter(ptr, name);
		if (boost::iequals(obj, "acquisitions"))
			return cGT_acquisitionsParameter(ptr, name);
		if (boost::iequals(obj, "gadget_chain")) {
			GadgetChain& gc = objectFromHandle<GadgetChain>(ptr);
			boost::shared_ptr<aGadget> sptr = gc.gadget_sptr(name);
			if (sptr.get())
				return sptrObjectHandle(sptr);
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
			boost::shared_ptr<CoilSensitivitiesContainer>
				csms(new CoilSensitivitiesAsImages(file));
			return sptrObjectHandle<CoilSensitivitiesContainer>(csms);
		}
		else {
			boost::shared_ptr<CoilSensitivitiesContainer>
				csms(new CoilSensitivitiesAsImages());
			return sptrObjectHandle<CoilSensitivitiesContainer>(csms);
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
		AcquisitionsContainer& acqs =
			objectFromHandle<AcquisitionsContainer>(h_acqs);
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
		AcquisitionsContainer& acqs =
			objectFromHandle<AcquisitionsContainer>(h_acqs);
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
		double* re = (double*)ptr_re;
		double* im = (double*)ptr_im;
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
	double* re = (double*)ptr_re;
	double* im = (double*)ptr_im;
	CAST_PTR(DataHandle, h_csms, ptr_csms);
	CoilDataContainer& list =
		objectFromHandle<CoilDataContainer>(h_csms);
	list.get_data(csm_num, re, im);
}

extern "C"
void
cGT_getCoilDataAbs(void* ptr_csms, int csm_num, size_t ptr)
{
	double* v = (double*)ptr;
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
		boost::shared_ptr<AcquisitionsContainer> acqs =
			objectSptrFromHandle<AcquisitionsContainer>(h_acqs);
		boost::shared_ptr<ImagesContainer> imgs =
			objectSptrFromHandle<ImagesContainer>(h_imgs);
		boost::shared_ptr<AcquisitionModel> am(new AcquisitionModel(acqs, imgs));
		return sptrObjectHandle<AcquisitionModel>(am);
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
		AcquisitionModel& am = objectFromHandle<AcquisitionModel>(h_am);
		boost::shared_ptr<CoilSensitivitiesContainer> sptr_csms =
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
		AcquisitionModel& am = objectFromHandle<AcquisitionModel>(h_am);
		ImagesContainer& imgs = objectFromHandle<ImagesContainer>(h_imgs);
		boost::shared_ptr<AcquisitionsContainer> sptr_acqs = am.fwd(imgs);
		return sptrObjectHandle<AcquisitionsContainer>(sptr_acqs);
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
		AcquisitionModel& am = objectFromHandle<AcquisitionModel>(h_am);
		AcquisitionsContainer& acqs =
			objectFromHandle<AcquisitionsContainer>(h_acqs);
		boost::shared_ptr<ImagesContainer> sptr_imgs = am.bwd(acqs);
		return sptrObjectHandle<ImagesContainer>(sptr_imgs);
	}
	CATCH;
}

extern "C"
void*
cGT_orderAcquisitions(void* ptr_acqs)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		AcquisitionsContainer& acqs =
			objectFromHandle<AcquisitionsContainer>(h_acqs);
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
		boost::shared_ptr<AcquisitionsContainer> 
			acquisitions(new AcquisitionsFile(file));
		return sptrObjectHandle<AcquisitionsContainer>(acquisitions);
	}
	CATCH;
}

extern "C"
void*
cGT_ISMRMRDAcquisitionsFile(const char* file)
{
	try {
		boost::shared_ptr<AcquisitionsContainer> 
			acquisitions(new AcquisitionsFile(file, true, true));
		return sptrObjectHandle<AcquisitionsContainer>(acquisitions);
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
		AcquisitionsContainer& input =
			objectFromHandle<AcquisitionsContainer>(h_input);
		proc.process(input);
		boost::shared_ptr<AcquisitionsContainer> sptr_ac = proc.get_output();
		return sptrObjectHandle<AcquisitionsContainer>(sptr_ac);
	}
	CATCH;
}

extern "C"
void*
cGT_acquisitionFromContainer(void* ptr_acqs, unsigned int acq_num)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		AcquisitionsContainer& acqs =
			objectFromHandle<AcquisitionsContainer>(h_acqs);
		boost::shared_ptr<ISMRMRD::Acquisition>
			sptr_acq(new ISMRMRD::Acquisition);
		acqs.get_acquisition(acq_num, *sptr_acq);
		return sptrObjectHandle<ISMRMRD::Acquisition>(sptr_acq);
	}
	CATCH;
}

extern "C"
void*
cGT_getAcquisitionsDimensions(void* ptr_acqs, size_t ptr_dim)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		AcquisitionsContainer& acqs =
			objectFromHandle<AcquisitionsContainer>(h_acqs);
		boost::shared_ptr<ISMRMRD::Acquisition>
			sptr_acq(new ISMRMRD::Acquisition);
		int not_reg = acqs.get_acquisitions_dimensions(ptr_dim);
		return dataHandle(not_reg);
		//return new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cGT_getAcquisitionsData
(void* ptr_acqs, unsigned int slice, size_t ptr_re, size_t ptr_im)
{
	try {
		double* re = (double*)ptr_re;
		double* im = (double*)ptr_im;
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		AcquisitionsContainer& acqs =
			objectFromHandle<AcquisitionsContainer>(h_acqs);
		int n = acqs.get_acquisitions_data(slice, re, im);
		return dataHandle(n);
		//return new DataHandle;
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
	if (boost::iequals(name, "number_of_samples"))
		return dataHandle((int)acq.number_of_samples());
	if (boost::iequals(name, "active_channels"))
		return dataHandle((int)acq.active_channels());
	if (boost::iequals(name, "trajectory_dimensions"))
		return dataHandle((int)acq.trajectory_dimensions());
	if (boost::iequals(name, "flags"))
		return dataHandle((int)acq.flags());
	if (boost::iequals(name, "idx_kspace_encode_step_1"))
		return dataHandle((int)acq.idx().kspace_encode_step_1);
	if (boost::iequals(name, "idx_repetition"))
		return dataHandle((int)acq.idx().repetition);
	if (boost::iequals(name, "idx_slice"))
		return dataHandle((int)acq.idx().slice);
	return parameterNotFound(name, __FILE__, __LINE__);
}

extern "C"
void*
cGT_acquisitionsParameter(void* ptr_acqs, const char* name)
{
	try {
		CAST_PTR(DataHandle, h_acqs, ptr_acqs);
		AcquisitionsContainer& acqs =
			objectFromHandle<AcquisitionsContainer>(h_acqs);
		if (boost::iequals(name, "undersampled"))
			return dataHandle((int)acqs.undersampled());
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
		AcquisitionsContainer& input = objectFromHandle<AcquisitionsContainer>(h_input);
		recon.process(input);
		boost::shared_ptr<ImagesContainer> sptr_img = recon.get_output();
		return sptrObjectHandle<ImagesContainer>(sptr_img);
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
		boost::shared_ptr<ImagesContainer> sptr_img = recon.get_output();
		return sptrObjectHandle<ImagesContainer>(sptr_img);
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
		ImagesContainer& input = objectFromHandle<ImagesContainer>(h_input);
		proc.process(input);
		boost::shared_ptr<ImagesContainer> sptr_img = proc.get_output();
		return sptrObjectHandle<ImagesContainer>(sptr_img);
	}
	CATCH;

}

extern "C"
void*
cGT_selectImages(void* ptr_input, unsigned int inc, unsigned int off)
{
	try {
		CAST_PTR(DataHandle, h_input, ptr_input);
		ImagesContainer& input = objectFromHandle<ImagesContainer>(h_input);
		boost::shared_ptr<ImagesContainer> sptr_img = input.clone(inc, off);
		return sptrObjectHandle<ImagesContainer>(sptr_img);
	}
	CATCH;
}

extern "C"
void*
cGT_imagesCopy(const void* ptr_imgs)
{
	try {
		CAST_PTR(DataHandle, h_imgs, ptr_imgs);
		ImagesContainer& imgs = 
			(ImagesContainer&)objectFromHandle<ImagesContainer>(h_imgs);
		boost::shared_ptr<ImagesContainer> clone = imgs.clone();
		return sptrObjectHandle<ImagesContainer>(clone);
	}
	CATCH;
}

extern "C"
void*
cGT_writeImages(void* ptr_imgs, const char* out_file, const char* out_group)
{
	try {
		CAST_PTR(DataHandle, h_imgs, ptr_imgs);
		ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
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
	ImagesContainer& images = objectFromHandle<ImagesContainer>(h_imgs);
	return sptrObjectHandle<ImageWrap>(images.sptr_image_wrap(img_num));
}

extern "C"
void
cGT_setImageToRealConversion(void* ptr_imgs, int type)
{
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& images = objectFromHandle<ImagesContainer>(h_imgs);
	images.set_image_to_real_conversion(type);
}

extern "C"
void
cGT_getImageDimensions(void* ptr_imgs, int img_num, size_t ptr_dim)
{
	int* dim = (int*)ptr_dim;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	list.get_image_dimensions(img_num, dim);
}

extern "C"
void
cGT_getImageDataAsDoubleArray(void* ptr_imgs, int img_num, size_t ptr_data)
{
	double* data = (double*)ptr_data;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	list.get_image_data_as_double_array(img_num, data);
}

extern "C"
void
cGT_getImageDataAsComplexArray(void* ptr_imgs, int img_num, size_t ptr_data)
{
	complex_float_t* data = (complex_float_t*)ptr_data;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	list.get_image_data_as_complex_array(img_num, data);
}

extern "C"
void
cGT_getImageDataAsCmplxArray
(void* ptr_imgs, int img_num, double* re, double* im)
{
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	list.get_image_data_as_cmplx_array(img_num, re, im);
}

extern "C"
void
cGT_getImagesDataAsDoubleArray(void* ptr_imgs, size_t ptr_data)
{
	double* data = (double*)ptr_data;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	list.get_images_data_as_double_array(data);
}

extern "C"
void
cGT_getImagesDataAsComplexArray(void* ptr_imgs, size_t ptr_re, size_t ptr_im)
{
	double* re = (double*)ptr_re;
	double* im = (double*)ptr_im;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	list.get_images_data_as_complex_array(re, im);
}

extern "C"
void*
cGT_imageTypes(const void* ptr_x)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		ImagesContainer& x = objectFromHandle<ImagesContainer>(h_x);
		int* result = (int*)malloc(sizeof(int));
		*result = x.types();
		DataHandle* handle = new DataHandle;
		handle->set(result, 0, GRAB);
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void*
cGT_imageDataType(const void* ptr_x, int im_num)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		ImagesContainer& x = objectFromHandle<ImagesContainer>(h_x);
		int* result = (int*)malloc(sizeof(int));
		*result = x.image_data_type(im_num);
		DataHandle* handle = new DataHandle;
		handle->set(result, 0, GRAB);
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void*
cGT_dataItems(const void* ptr_x)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		aDataContainer& x = objectFromHandle<aDataContainer>(h_x);
		int* result = (int*)malloc(sizeof(int));
		*result = x.items();
		DataHandle* handle = new DataHandle;
		handle->set(result, 0, GRAB);
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void*
cGT_norm(const void* ptr_x)
{
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		aDataContainer& x = objectFromHandle<aDataContainer>(h_x);
		double* result = (double*)malloc(sizeof(double));
		*result = x.norm();
		DataHandle* handle = new DataHandle;
		handle->set(result, 0, GRAB);
		return (void*)handle;
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
		aDataContainer& x = objectFromHandle<aDataContainer>(h_x);
		aDataContainer& y = objectFromHandle<aDataContainer>(h_y);
		complex_double_t* result =
			(complex_double_t*)malloc(sizeof(complex_double_t));
		*result = x.dot(y);
		DataHandle* handle = new DataHandle;
		handle->set(result, 0, GRAB);
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void*
cGT_axpby(
double ar, double ai, const void* ptr_x,
double br, double bi, const void* ptr_y
){
	try {
		CAST_PTR(DataHandle, h_x, ptr_x);
		CAST_PTR(DataHandle, h_y, ptr_y);
		aDataContainer& x = objectFromHandle<aDataContainer>(h_x);
		aDataContainer& y = objectFromHandle<aDataContainer>(h_y);
		boost::shared_ptr<aDataContainer> sptr_z = x.new_data_container();
		complex_double_t a(ar, ai);
		complex_double_t b(br, bi);
		sptr_z->axpby(a, x, b, y);
		return sptrObjectHandle<aDataContainer>(sptr_z);
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
		boost::shared_ptr<aGadget>& g = objectSptrFromHandle<aGadget>(h_r);
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
		boost::shared_ptr<aGadget>& g = objectSptrFromHandle<aGadget>(h_w);
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
		boost::shared_ptr<aGadget>& g = objectSptrFromHandle<aGadget>(h_g);
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
		size_t i, j, k;
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
		boost::shared_ptr<ImagesContainer> sptr_images =
			objectSptrFromHandle<ImagesContainer>(h_img);
		con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
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
		ImagesContainer& images = objectFromHandle<ImagesContainer>(h_img);
		for (int i = 0; i < images.number(); i++) {
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

