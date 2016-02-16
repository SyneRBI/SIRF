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
#include "xgadgetron.h"
#include "gadget_lib.h"

boost::shared_ptr<boost::mutex> Mutex::sptr_mutex_;
//boost::shared_ptr<GadgetronClientConnector> GTConnector::sptr_con_;

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
		//else if (boost::iequals(name, "MRIReconstruction"))
		//	return newObjectHandle<GadgetChain, MRIReconstruction>();
		else if (boost::iequals(name, "ImageReconstructor"))
			return newObjectHandle<GadgetChain, ImageReconstructor>();
		else if (boost::iequals(name, "ImagesProcessor"))
			return newObjectHandle<GadgetChain, ImagesProcessor>();
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
	CATCH
}

extern "C"
void*
cGT_ISMRMRDAcquisitionsFromFile(const char* file)
{
	try {
		boost::shared_ptr<AcquisitionsContainer> acquisitions(new AcquisitionsFile(file));
		return sptrObjectHandle<AcquisitionsContainer>(acquisitions);
	}
	CATCH
}

extern "C"
void*
cGT_acquisitionsProcessor(const char* file)
{
	try {
		boost::shared_ptr<AcquisitionsProcessor> proc(new AcquisitionsProcessor(file));
		return sptrObjectHandle<AcquisitionsProcessor>(proc);
	}
	CATCH
}

extern "C"
void*
cGT_ISMRMRDatasetFromFile(const char* file, const char* group)
{
	try {
		boost::shared_ptr<ISMRMRD::Dataset> ismrmrd_dataset
			(new ISMRMRD::Dataset(file, group, false));
		return sptrObjectHandle<ISMRMRD::Dataset>(ismrmrd_dataset);
	}
	CATCH
}

extern "C"
void*
cGT_readISMRMRDatasetHeader(void* ptr_data, void* ptr_head)
{
	try {
		CAST_PTR(DataHandle, h_data, ptr_data);
		CAST_PTR(DataHandle, h_head, ptr_head);
		ISMRMRD::Dataset& data = objectFromHandle<ISMRMRD::Dataset>(h_data);
		std::string& head = objectFromHandle<std::string>(h_head);
		data.readHeader(head);
	}
	CATCH

	return (void*)new DataHandle;
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
	CATCH

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
	CATCH

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
	CATCH

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
	CATCH

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
	CATCH

		return (void*)new DataHandle;
}

extern "C"
void*
cGT_setEndGadget(void* ptr_gc, const void* ptr_g)
{
	try {
		CAST_PTR(DataHandle, h_gc, ptr_gc);
		CAST_PTR(DataHandle, h_g, ptr_g);
		GadgetChain& gc = objectFromHandle<GadgetChain>(h_gc);
		boost::shared_ptr<aGadget>& g = objectSptrFromHandle<aGadget>(h_g);
		gc.set_endgadget(g);
	}
	CATCH

		return (void*)new DataHandle;
}

extern "C"
void*
cGT_reconstructImages(void* ptr_recon, void* ptr_input)
{
	try {
		CAST_PTR(DataHandle, h_recon, ptr_recon);
		CAST_PTR(DataHandle, h_input, ptr_input);
		//MRIReconstruction& recon = objectFromHandle<MRIReconstruction>(h_recon);
		ImageReconstructor& recon = objectFromHandle<ImageReconstructor>(h_recon);
		AcquisitionsContainer& input = objectFromHandle<AcquisitionsContainer>(h_input);
		recon.process(input);
		boost::shared_ptr<ImagesContainer> sptr_im = recon.get_output();
		ObjectHandle<ImagesContainer>* ptr_handle = new ObjectHandle<ImagesContainer>(sptr_im);
		return (void*)ptr_handle;
	}
	CATCH;

}

extern "C"
void*
cGT_reconstructedImagesList(void* ptr_recon)
{
	try {
		CAST_PTR(DataHandle, h_recon, ptr_recon);
		//MRIReconstruction& recon = objectFromHandle<MRIReconstruction>(h_recon);
		ImageReconstructor& recon = objectFromHandle<ImageReconstructor>(h_recon);
		boost::shared_ptr<ImagesContainer> sptr_im = recon.get_output();
		ObjectHandle<ImagesContainer>* ptr_handle = new ObjectHandle<ImagesContainer>(sptr_im);
		return (void*)ptr_handle;
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
		boost::shared_ptr<ImagesContainer> sptr_im = proc.get_output();
		ObjectHandle<ImagesContainer>* ptr_handle = 
			new ObjectHandle<ImagesContainer>(sptr_im);
		return (void*)ptr_handle;
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
		boost::shared_ptr<AcquisitionsContainer> sptr_im = proc.get_output();
		ObjectHandle<AcquisitionsContainer>* ptr_handle = 
			new ObjectHandle<AcquisitionsContainer>(sptr_im);
		return (void*)ptr_handle;
	}
	CATCH;

}

//extern "C"
//void*
//cGT_processedImagesList(void* ptr_proc)
//{
//	try {
//		CAST_PTR(DataHandle, h_proc, ptr_proc);
//		ImageProcessor& proc = objectFromHandle<ImageProcessor>(h_proc);
//		boost::shared_ptr<ImagesList> sptr_im = proc.get_output();
//		ObjectHandle<ImagesList>* ptr_handle = new ObjectHandle<ImagesList>(sptr_im);
//		return (void*)ptr_handle;
//	}
//	CATCH;
//
//}

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
	CATCH

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
	CATCH

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
	CATCH

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
	CATCH

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
	CATCH

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
	CATCH

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
		//boost::mutex& mtx = conn.mutex();
		con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE,
			boost::shared_ptr<GadgetronClientMessageReader>
			(new GadgetronClientImageMessageReader(file, group, &mtx)));
	}
	CATCH

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
	CATCH

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
		//boost::mutex& mtx = conn.mutex();
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
	CATCH

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
			ImageWrap& iw = images.imageWrap(i);
			con.send_wrapped_image(iw);
		}
	}
	CATCH

		return (void*)new DataHandle;
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
	CATCH

	return (void*)new DataHandle;
}

extern "C"
int
cGT_numImages(void* ptr_imgs)
{
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	return list.number();
}

extern "C"
void
cGT_getImageDimensions(void* ptr_imgs, int im_num, size_t ptr_dim)
{
	int* dim = (int*)ptr_dim;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	list.getImageDimensions(im_num, dim);
}

extern "C"
void
cGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, size_t ptr_data)
{
	double* data = (double*)ptr_data;
	CAST_PTR(DataHandle, h_imgs, ptr_imgs);
	ImagesContainer& list = objectFromHandle<ImagesContainer>(h_imgs);
	list.getImageDataAsDoubleArray(im_num, data);
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
	CATCH

	return (void*)new DataHandle;
}