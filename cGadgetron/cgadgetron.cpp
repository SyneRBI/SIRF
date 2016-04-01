#include <boost/program_options.hpp>
#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>

#include <fstream>
#include <streambuf>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <exception>
#include <map>
#include <thread>
#include <chrono>
#include <condition_variable>

#include "gadgetron_data_containers.h"
#include "gadgetron_client.h"
#include "data_handle.h"
#include "iutilities.h"
#include "cgadgetron_par.h"
#include "xgadgetron.h"
#include "gadget_lib.h"
#include "chain_lib.h"

#define GRAB 1

boost::shared_ptr<boost::mutex> Mutex::sptr_mutex_;

static void*
unknownObject(const char* obj, const char* name, const char* file, int line)
{
	DataHandle* handle = new DataHandle;
	std::string error = "unknown ";
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
	std::string error = "parameter ";
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
		else if (boost::iequals(name, "GTConnector"))
			return newObjectHandle<GTConnector, GTConnector>();
		else if (boost::iequals(name, "string"))
			return newObjectHandle<std::string, std::string>();
		else if (boost::iequals(name, "ImagesList"))
			return newObjectHandle<ImagesContainer, ImagesList>();
		else if (boost::iequals(name, "GadgetChain"))
			return newObjectHandle<GadgetChain, GadgetChain>();
		else if (boost::iequals(name, "AcquisitionsProcessor"))
			return newObjectHandle<GadgetChain, AcquisitionsProcessor>();
		else if (boost::iequals(name, "ImagesReconstructor"))
			return newObjectHandle<GadgetChain, ImagesReconstructor>();
		else if (boost::iequals(name, "ImagesProcessor"))
			return newObjectHandle<GadgetChain, ImagesProcessor>();
		else if (boost::iequals(name, "RemoveOversamplingProcessor"))
			return newObjectHandle<GadgetChain, RemoveOversamplingProcessor>();
		else if (boost::iequals(name, "SimpleReconstructionProcessor"))
			return newObjectHandle<GadgetChain, SimpleReconstructionProcessor>();
		else if (boost::iequals(name, "ExtractRealImagesProcessor"))
			return newObjectHandle<GadgetChain, ExtractRealImagesProcessor>();
		else if (boost::iequals(name, "GadgetIsmrmrdAcquisitionMessageReader"))
			return newObjectHandle<aGadget, IsmrmrdAcqMsgReader>();
		else if (boost::iequals(name, "GadgetIsmrmrdAcquisitionMessageWriter"))
			return newObjectHandle<aGadget, IsmrmrdAcqMsgWriter>();
		else if (boost::iequals(name, "MRIImageReader"))
			return newObjectHandle<aGadget, IsmrmrdImgMsgReader>();
		else if (boost::iequals(name, "MRIImageWriter"))
			return newObjectHandle<aGadget, IsmrmrdImgMsgWriter>();
		else if (boost::iequals(name, "RemoveROOversamplingGadget"))
			return newObjectHandle<aGadget, RemoveOversamplingGadget>();
		else if (boost::iequals(name, "AcquisitionAccumulateTriggerGadget"))
			return newObjectHandle<aGadget, AcqAccTrigGadget>();
		else if (boost::iequals(name, "BucketToBufferGadget"))
			return newObjectHandle<aGadget, BucketToBuffGadget>();
		else if (boost::iequals(name, "SimpleReconGadget"))
			return newObjectHandle<aGadget, SimpleReconstructionGadget>();
		else if (boost::iequals(name, "ImageArraySplitGadget"))
			return newObjectHandle<aGadget, ImgArrSplitGadget>();
		else if (boost::iequals(name, "ExtractGadget"))
			return newObjectHandle<aGadget, ExtGadget>();
		else if (boost::iequals(name, "ImageFinishGadget"))
			return newObjectHandle<aGadget, ImgFinishGadget>();
		else if (boost::iequals(name, "AcquisitionFinishGadget"))
			return newObjectHandle<aGadget, AcqFinishGadget>();
		else if (boost::iequals(name, "SimpleReconGadgetSet"))
			return newObjectHandle<aGadget, SimpleReconstructionGadgetSet>();
		std::cout << "object " << name << "		not found" << std::endl;
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
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void*
cGT_CoilSensitivitiesFromFile(const char* file)
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
void
cGT_getCSMDimensions(void* ptr_csms, int csm_num, size_t ptr_dim)
{
	int* dim = (int*)ptr_dim;
	CAST_PTR(DataHandle, h_csms, ptr_csms);
	CoilSensitivitiesContainer& list = 
		objectFromHandle<CoilSensitivitiesContainer>(h_csms);
	list.get_dim(csm_num, dim);
}

extern "C"
void
cGT_getCSMData(void* ptr_csms, int csm_num, size_t ptr_re, size_t ptr_im)
{
	double* re = (double*)ptr_re;
	double* im = (double*)ptr_im;
	CAST_PTR(DataHandle, h_csms, ptr_csms);
	CoilSensitivitiesContainer& list =
		objectFromHandle<CoilSensitivitiesContainer>(h_csms);
	list.get_data(csm_num, re, im);
}

extern "C"
void
cGT_getCSMDataAbs(void* ptr_csms, int csm_num, size_t ptr)
{
	double* v = (double*)ptr;
	CAST_PTR(DataHandle, h_csms, ptr_csms);
	CoilSensitivitiesContainer& list =
		objectFromHandle<CoilSensitivitiesContainer>(h_csms);
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
	else
		return parameterNotFound(name, __FILE__, __LINE__);
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
		if (boost::iequals(g.name(), "AcquisitionAccumulateTriggerGadget"))
			g.set_property(prop, value);
		else if (boost::iequals(g.name(), "BucketToBufferGadget"))
			g.set_property(prop, value);
		else if (boost::iequals(g.name(), "SimpleReconGadgetSet"))
			g.set_property(prop, value);
		//else
		//	return unknownObject
		//		("gadget with properties", g.name().c_str(), __FILE__, __LINE__);
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
cGT_registerHDFReceiver(void* ptr_con, const char* file, const char* group)
{
	try {
		CAST_PTR(DataHandle, h_con, ptr_con);
		GTConnector& conn = objectFromHandle<GTConnector>(h_con);
		GadgetronClientConnector& con = conn();
		Mutex mutex;
		boost::mutex& mtx = mutex();
		con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageReader(file, group, &mtx)));
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

