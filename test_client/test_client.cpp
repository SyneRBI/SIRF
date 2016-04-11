/*****************************************
*  Standalone ISMRMRD Gadgetron Client
*
* Author: Michael S. Hansen
*
* Dependencies: ISMRMRD and Boost
*
*****************************************/

//TODO:
// -Blobs (for DICOM image support)
//  - First implementation is in, but testing needed
// -NIFTI and Analyze output
// -Check on potential threading problem with asio socket 
//    - having and reading and writing thread is supposedly not safe, but seems to work here
// -Add command line switch for controlling verbosity of output
// -Static linking for standalone executable. 

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

#include "gadgetron_client.h"
#include "data_handle.h"
#include "cgadgetron.h"
#include "xgadgetron.h"
#include "gadget_lib.h"
#include "gadgetron_data_containers.h"
#include "iutilities.h"

#include "text_writer.h"

//aTextWriter* TextWriterHandle::information_channel_;
//aTextWriter* TextWriterHandle::warning_channel_;
//aTextWriter* TextWriterHandle::error_channel_;

std::string get_date_time_string()
{
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	//std::cout << rawtime << std::endl;
	//std::cout << clock() << std::endl;

//	auto now = std::chrono::system_clock::now();
//	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
////		- std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch());
//	std::cout << ms.count() << " ms" << std::endl;
	std::cout << xGadgetronUtilities::milliseconds() << " ms" << std::endl;

	//struct timeval tp;
	//gettimeofday(&tp, NULL);
	//long long mslong = (long long)tp.tv_sec * 1000L + tp.tv_usec / 1000;
	//std::cout << mslong << std::endl;

	std::stringstream str;
	str << timeinfo->tm_year + 1900 << "-"
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_mon + 1 << "-"
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_mday << " "
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_hour << ":"
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_min << ":"
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_sec;

	std::string ret = str.str();

	return ret;
}

//int test5(
//	const char* host,
//	const char* port,
//	unsigned int timeout,
//	const char* in_file,
//	const char* in_group,
//	const char* out_file,
//	const char* out_group);

int test6(
	const char* host,
	const char* port,
	unsigned int timeout,
	const char* in_file,
	const char* in_group,
	const char* out_file,
	const char* out_group);

int test7(
	const char* host,
	const char* port,
	unsigned int timeout,
	const char* in_file,
	const char* in_group,
	const char* out_file,
	const char* out_group);

#define NTRIALS 5

int test8(const char* in_file);

//int test9(const char* in_file);

int test10();

int test11(const char* in_file);

int test12();

int test13(const char* in_file);

namespace po = boost::program_options;
using boost::asio::ip::tcp;

int main(int argc, char **argv)
{
	std::string host_name;
	std::string port;
	std::string in_filename;
	std::string out_filename;
	std::string hdf5_in_group;
	std::string hdf5_out_group;
	std::string config_file;
	std::string config_file_local;
	std::string config_xml_local;
	unsigned int loops;
	unsigned int timeout_ms;
	std::string out_fileformat;
	bool open_input_file = true;

	void* printer = newTextPrinter("stdout");
	//void* writer = newTextWriter("stdout.txt");
	openChannel(0, printer);

	po::options_description desc("Allowed options");

	desc.add_options()
		("help,h", "produce help message")
		("query,q", "Dependency query mode")
		("port,p", po::value<std::string>(&port)->default_value("9002"), "Port")
		("address,a", po::value<std::string>(&host_name)->default_value("localhost"), "Address (hostname) of Gadgetron host")
		("filename,f", po::value<std::string>(&in_filename), "Input file")
		("outfile,o", po::value<std::string>(&out_filename)->default_value("out.h5"), "Output file")
		("in-group,g", po::value<std::string>(&hdf5_in_group)->default_value("/dataset"), "Input data group")
		("out-group,G", po::value<std::string>(&hdf5_out_group)->default_value(get_date_time_string()), "Output group name")
		("config,c", po::value<std::string>(&config_file)->default_value("default.xml"), "Configuration file (remote)")
		("config-local,C", po::value<std::string>(&config_file_local), "Configuration file (local)")
		("loops,l", po::value<unsigned int>(&loops)->default_value(1), "Loops")
		("timeout,t", po::value<unsigned int>(&timeout_ms)->default_value(10000), "Timeout [ms]");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	}

	if (!vm.count("filename") && !vm.count("query")) {
		std::cout << std::endl << std::endl << "\tYou must supply a filename" << std::endl << std::endl;
		std::cout << desc << std::endl;
		return -1;
	}

	if (vm.count("query")) {
		open_input_file = false;
	}

	//return test13(in_filename.c_str());

	//return test12();

	//return test11(in_filename.c_str());

	return test10();

	//return test9(in_filename.c_str());

	//return test8(in_filename.c_str());

	return test7(
		host_name.c_str(),
		port.c_str(),
		timeout_ms,
		in_filename.c_str(),
		hdf5_in_group.c_str(),
		"out7.h5",
		hdf5_out_group.c_str());

	//return test6(
	//	host_name.c_str(),
	//	port.c_str(),
	//	timeout_ms,
	//	in_filename.c_str(),
	//	hdf5_in_group.c_str(),
	//	"out6.h5",
	//	hdf5_out_group.c_str());

	//return test5(
	//	host_name.c_str(),
	//	port.c_str(),
	//	timeout_ms,
	//	in_filename.c_str(),
	//	hdf5_in_group.c_str(),
	//	"out5.h5",
	//	hdf5_out_group.c_str());

	return 0;
}

#if 0
int test5(
	const char* host,
	const char* port,
	unsigned int timeout,
	const char* in_file,
	const char* in_group,
	const char* out_file,
	const char* out_group)
{
	std::cout << "Gadgetron ISMRMRD client" << std::endl;
	std::cout << "  -- host            :      " << host << std::endl;
	std::cout << "  -- port            :      " << port << std::endl;
	std::cout << "  -- hdf5 file  in   :      " << in_file << std::endl;
	std::cout << "  -- hdf5 group in   :      " << in_group << std::endl;
	std::cout << "  -- hdf5 file out   :      " << out_file << std::endl;
	std::cout << "  -- hdf5 group out  :      " << out_group << std::endl;

	try {
		//ISMRMRD::Dataset input(in_file, in_group, false);
		boost::shared_ptr<AcquisitionsContainer> sptr_input(new AcquisitionsFile(in_file));
		AcquisitionsContainer& input = *sptr_input;
		//AcquisitionsFile input(in_file);
		//std::cout << "ok" << std::endl;

		SPTR(aGadget, ro, RemoveOversamplingGadget);
		SPTR(aGadget, aat, AcqAccTrigGadget);
		SPTR(aGadget, bb, BucketToBuffGadget);
		SPTR(aGadget, sr, SimpleReconstructionGadget);
		SPTR(aGadget, ias, ImgArrSplitGadget);
		SPTR(aGadget, ext, ExtGadget);
		SPTR(aGadget, fin, ImgFinishGadget);
		//std::cout << "ok" << std::endl;

		{
			ImagesReconstructor recon;
			//std::cout << "ok" << std::endl;

			recon.add_gadget("g1", ro);
			recon.add_gadget("g2", aat);
			recon.add_gadget("g3", bb);
			recon.add_gadget("g4", sr);
			recon.add_gadget("g5", ias);
			//recon.add_gadget("g6", ext);
			//recon.add_gadget("g7", fin);
			//std::cout << "ok" << std::endl;

			recon.process(input);
			//std::cout << "ok" << std::endl;

			ImagesList& imgs = (ImagesList&)*recon.get_output();
			//std::cout << "ok" << std::endl;

			ImageWrap& iw_cmplx = imgs.image_wrap(0);
			//input.getPhantomAsComplexFloat(iw_cmplx);
			//std::cout << iw_cmplx.dot(iw_cmplx) << std::endl;

			//ImageHandle& ih_cmplx = imgs.image_handle(0);
			//std::cout << ih_cmplx.dot(ih_cmplx) << std::endl;

			ImageWrap& iw_cmplx1 = imgs.image_wrap(1);
			//ImageHandle& ih_cmplx1 = imgs.image_handle(1);
			//input.getPhantomAsComplexFloat(iw_cmplx1);

			std::cout << imgs.dot(imgs) << std::endl;

			ImageWrap iw_cmplx2(iw_cmplx1);
			complex_float_t a = -1.0;
			complex_float_t b = 1.0;
			iw_cmplx2.axpby(a, iw_cmplx, b);
			//std::cout << iw_cmplx2.norm() << std::endl;
			//std::cout << iw_cmplx.diff(iw_cmplx1) << std::endl;

			ImagesList list(imgs);
			//list.axpby(a, imgs, b);
			//std::cout << list.norm() << std::endl;

			ImagesProcessor proc;
			proc.add_gadget("g1", ext);
			//proc.add_gadget("g2", fin);
			proc.process(imgs);
			boost::shared_ptr<ImagesContainer> sptr_imgs = proc.get_output();
			ImagesList& images = (ImagesList&)*sptr_imgs;
			//ImagesList& images = (ImagesList&)*proc.get_output();

			//std::cout << images.dot(images) << std::endl;
			//std::cout << images.number() << std::endl;

			//images.write(out_file, out_group, GTConnector());
			images.write(out_file, out_group);
			//std::cout << "ok" << std::endl;

			//AcquisitionModel acq_mod(input);
			AcquisitionModel acq_mod(sptr_input, sptr_imgs);

			ImageWrap& iw = images.image_wrap(0);

			//input.getPhantomAsFloat(iw);
			//std::cout << iw.dot(iw) << std::endl;
			std::cout << images.dot(images) << std::endl;

			//images.write(out_file, out_group);

			std::string acq_file = out_file;
			acq_file += out_group;
			boost::replace_all(acq_file, " ", "_");
			boost::replace_all(acq_file, ":", "_");
			boost::replace_all(acq_file, ".h5", "_");
			acq_file += ".h5";
			std::cout << acq_file << std::endl;

			AcquisitionsFile acqs(acq_file, true, true);

			//ImageHandle& ih = images.image_handle(0);

			acq_mod.fwd(imgs, acqs);
			//acq_mod.fwd(iw, acqs);
			//acq_mod.fwd(ih, acqs);

			//ISMRMRD::Acquisition a;
			//ISMRMRD::Acquisition b;
			//input.getAcquisition(0, a);
			//acqs.getAcquisition(0, b);
			//std::cout << AcquisitionsContainer::diff(a, b) << std::endl;

			std::cout << input.diff(acqs) << std::endl;

			//std::cout << input.dot(input) << std::endl;
			std::cout << acqs.dot(input) << std::endl;
			std::cout << acqs.dot(acqs) << std::endl;
			std::cout << pow(acqs.norm(), 2) << std::endl;

			//ISMRMRD::Acquisition acq;
			//acqs.getAcquisition(0, acq);
			//std::cout << AcquisitionsContainer::dot(acq, acq) << std::endl;

			//acq_mod.bwd(ih_cmplx, acqs);
			//std::cout << ih_cmplx.dot(ih_cmplx1) << std::endl;
			//acq_mod.bwd(iw_cmplx, acqs);
			//std::cout << iw_cmplx.dot(iw_cmplx1) << std::endl;
			ImagesList il(imgs);
			acq_mod.bwd(imgs, acqs);
			std::cout << imgs.dot(il) << std::endl;

			//AcquisitionsFile diff("tmp.h5", true, true);
			boost::shared_ptr<AcquisitionsContainer> sptr_diff =
				input.new_acquisitions_container();
			AcquisitionsContainer& diff = *sptr_diff;
			a = -acqs.dot(input) / input.dot(input);
			//AcquisitionsContainer::axpby(a, input, b, acqs, diff);
			diff.axpby(a, input, b, acqs);
			std::cout << diff.norm()/acqs.norm() << std::endl;
			//std::remove("tmp.h5");
		}
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}
#endif

int test6(
	const char* host,
	const char* port,
	unsigned int timeout,
	const char* in_file,
	const char* in_group,
	const char* out_file,
	const char* out_group)
{
	std::cout << "Gadgetron ISMRMRD client" << std::endl;
	std::cout << "  -- host            :      " << host << std::endl;
	std::cout << "  -- port            :      " << port << std::endl;
	std::cout << "  -- hdf5 file  in   :      " << in_file << std::endl;
	std::cout << "  -- hdf5 group in   :      " << in_group << std::endl;
	std::cout << "  -- hdf5 file out   :      " << out_file << std::endl;
	std::cout << "  -- hdf5 group out  :      " << out_group << std::endl;

	try {
		//void* h_data = cGT_ISMRMRDatasetFromFile(in_file, in_group);
		void* h_data = cGT_ISMRMRDAcquisitionsFromFile(in_file);

		void* h_ro = cGT_newObject("RemoveROOversamplingGadget");
		void* h_aat = cGT_newObject("AcquisitionAccumulateTriggerGadget");
		void* h_bb = cGT_newObject("BucketToBufferGadget");
		void* h_sr = cGT_newObject("SimpleReconGadget");
		void* h_ias = cGT_newObject("ImageArraySplitGadget");
		void* h_e = cGT_newObject("ExtractGadget");
		//void* h_if = cGT_newObject("ImageFinishGadget");

		void* h_recon = cGT_newObject("ImageReconstructor");
		cGT_addGadget(h_recon, "g1", h_ro);
		cGT_addGadget(h_recon, "g2", h_aat);
		cGT_addGadget(h_recon, "g3", h_bb);
		cGT_addGadget(h_recon, "g4", h_sr);
		cGT_addGadget(h_recon, "g5", h_ias);
		//cGT_addGadget(h_recon, "g6", h_e);
		//cGT_addGadget(h_recon, "g7", h_if);

		//void* h_images = cGT_runMRIReconstruction(h_recon, h_data);
		void* h_images = cGT_reconstructImages(h_recon, h_data);

		//void* h_images = cGT_reconstructedImagesList(h_recon);

		void* h_proc = cGT_newObject("ImagesProcessor");
		cGT_addGadget(h_proc, "g1", h_e);
		//cGT_addGadget(h_proc, "g2", h_if);
		void* h_imgs = cGT_processImages(h_proc, h_images);
		if (executionStatus(h_imgs)) {
			std::cout << "exception thrown" << std::endl;
			exit(1);
		}

		//std::cout << cGT_numImages(h_imgs) << std::endl;
		std::cout << cGT_dataItems(h_imgs) << std::endl;

		cGT_writeImages(h_imgs, out_file, out_group);

		deleteObject(h_ro);
		deleteObject(h_aat);
		deleteObject(h_bb);
		deleteObject(h_sr);
		deleteObject(h_ias);
		deleteObject(h_e);
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}

int test7(
	const char* host,
	const char* port,
	unsigned int timeout,
	const char* in_file,
	const char* in_group,
	const char* out_file,
	const char* out_group)
{
	std::cout << "Gadgetron ISMRMRD client" << std::endl;
	std::cout << "  -- host            :      " << host << std::endl;
	std::cout << "  -- port            :      " << port << std::endl;
	std::cout << "  -- hdf5 file  in   :      " << in_file << std::endl;
	std::cout << "  -- hdf5 group in   :      " << in_group << std::endl;
	std::cout << "  -- hdf5 file out   :      " << out_file << std::endl;
	std::cout << "  -- hdf5 group out  :      " << out_group << std::endl;

	std::cout << xGadgetronUtilities::milliseconds() << std::endl;

	try {
		void* h_input = cGT_ISMRMRDAcquisitionsFromFile(in_file);

		//DataHandle* h = (DataHandle*)h_input;
		//printf("%p\n", h->data());
		//std::cout << cGT_norm(h_input) << std::endl;

		void* h_ro = cGT_newObject("RemoveROOversamplingGadget");
		void* h_sr = cGT_newObject("SimpleReconGadgetSet");
		void* h_e = cGT_newObject("ExtractGadget");

		std::string acq_file = out_file;
		acq_file += out_group;
		boost::replace_all(acq_file, " ", "_");
		boost::replace_all(acq_file, ":", "_");
		boost::replace_all(acq_file, ".h5", "_");
		acq_file += ".h5";
		std::cout << acq_file << std::endl;
		void* h_proc = cGT_newObject("AcquisitionsProcessor");//cGT_acquisitionsProcessor();
		//void* h_proc = cGT_acquisitionsProcessor(acq_file.c_str());
		cGT_addGadget(h_proc, "g1", h_ro);

		void* h_output = cGT_processAcquisitions(h_proc, h_input);
		if (executionStatus(h_output)) {
			std::cout << "exception thrown" << std::endl;
			exit(1);
		}
		std::cout << acq_file << std::endl;
		void* h_output1 = cGT_processAcquisitions(h_proc, h_input);
		if (executionStatus(h_output)) {
			std::cout << "exception thrown" << std::endl;
			exit(1);
		}

		void* h_recon = cGT_newObject("ImagesReconstructor");
		cGT_addGadget(h_recon, "g1", h_sr);
		
		//cGT_addGadget(h_recon, "g2", h_e);

		void* h_images = cGT_reconstructImages(h_recon, h_output);

		void* h_proc_img = cGT_newObject("ImagesProcessor");
		cGT_addGadget(h_proc_img, "g1", h_e);
		void* h_imgs = cGT_processImages(h_proc_img, h_images);
		if (executionStatus(h_imgs)) {
			std::cout << "exception thrown" << std::endl;
			exit(1);
		}

		cGT_writeImages(h_imgs, out_file, out_group);

		deleteObject(h_ro);
		deleteObject(h_sr);
		deleteObject(h_e);
		deleteObject(h_images);
		deleteObject(h_imgs);
		deleteObject(h_recon);
		deleteObject(h_input);
		deleteObject(h_output);
		deleteObject(h_output1);
		deleteObject(h_proc);
		deleteObject(h_proc_img);
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}

int test8(const char* in_file)
{
	//GTConnector conn;
	//std::string filename("tmp.h5");
	//boost::shared_ptr<AcquisitionsContainer> sptr_acqs_;
	//sptr_acqs_.reset(new AcquisitionsFile(filename, true, true));
	//conn().register_reader(GADGET_MESSAGE_ISMRMRD_ACQUISITION,
	//	boost::shared_ptr<GadgetronClientMessageReader>
	//	(new GadgetronClientAcquisitionMessageCollector(sptr_acqs_)));

	boost::shared_ptr<AcquisitionsContainer> sptr_input(new AcquisitionsFile(in_file));
	AcquisitionsContainer& input = *sptr_input;

	//AcquisitionsProcessor proc;
	//proc.process(input);
	//boost::shared_ptr<AcquisitionsContainer> sptr_ac = proc.get_output();

	//void* ptr_input = cGT_ISMRMRDAcquisitionsFromFile(in_file);
	void* ptr_proc = cGT_newObject("AcquisitionsProcessor");
	//void* ptr_proc = cGT_acquisitionsProcessor();
	//void* h_acqs = cGT_processAcquisitions(h_proc, h_input);
	CAST_PTR(DataHandle, h_proc, ptr_proc);
	//CAST_PTR(DataHandle, h_input, ptr_input);
	AcquisitionsProcessor& proc =
		objectFromHandle<AcquisitionsProcessor>(h_proc);
	//AcquisitionsContainer& input =
	//	objectFromHandle<AcquisitionsContainer>(h_input);
	proc.process(input);
	//boost::shared_ptr<AcquisitionsContainer> sptr_ac = proc.get_output();
	////return sptrObjectHandle<AcquisitionsContainer>(sptr_ac);
	//ObjectHandle<AcquisitionsContainer>* ptr_handle =
	//	new ObjectHandle<AcquisitionsContainer>(sptr_ac);

	//deleteObject(ptr_input);
	deleteDataHandle(ptr_proc);
	//deleteObject(ptr_proc);
	//deleteObject(h_acqs);
	return 0;
}

#if 0
int test9(const char* in_file)
{
	std::cout << "Gadgetron ISMRMRD client" << std::endl;
	std::cout << "  -- hdf5 file  in   :      " << in_file << std::endl;

	try {
		boost::shared_ptr<AcquisitionsContainer> sptr_input(new AcquisitionsFile(in_file));
		AcquisitionsContainer& input = *sptr_input;

		//aVector<complex_double_t>& input_vector = (aVector<complex_double_t>&)input;
		//std::cout << input_vector.norm() << std::endl;

		//SPTR(aGadget, ro, RemoveOversamplingGadget);
		//SPTR(aGadget, sr, SimpleReconstructionGadgetSet);

		//ImagesReconstructor recon;
		//recon.add_gadget("g1", ro);
		//recon.add_gadget("g2", sr);
		//recon.process(input);

		//ImagesList& imgs = (ImagesList&)*recon.get_output();
		//ImageWrap iw(imgs.image_wrap(0));
		//ISMRMRD::Image<complex_float_t>* ptr_img = 
		//	(ISMRMRD::Image<complex_float_t>*)iw.ptr_image();
		//ISMRMRD::Image<complex_float_t>& img = *ptr_img;

		boost::shared_ptr<ISMRMRD::NDArray<complex_float_t> > sptr_csm = input.coils();
		int ndim = sptr_csm->getNDim();
		const size_t *dims = sptr_csm->getDims();
		for (int i = 0; i < ndim; i++)
			std::cout << dims[i] << ' ';
		std::cout << '\n';

		int nx = dims[0];
		int ny = dims[1];
		int nc = dims[2];

		//ISMRMRD::Image<complex_float_t>* ptr_img = 
		//	new ISMRMRD::Image<complex_float_t>(nx, ny, 1, nc);
		//ISMRMRD::Image<complex_float_t>& img = *ptr_img;
		CFImage* ptr_img = new CFImage(nx, ny, 1, nc);
		CFImage& img = *ptr_img;

		std::cout << img.getImageSeriesIndex() << std::endl;

		for (int ic = 0; ic < nc; ic++)
			for (int iy = 0; iy < ny; iy++)
				for (int ix = 0; ix < nx; ix++)
					img(ix, iy, 0, ic) = (*sptr_csm)(ix, iy, ic);

		//ImageWrap iw(ISMRMRD::ISMRMRD_CXFLOAT, ptr_img);
		Mutex mtx;
		mtx.lock();
		ISMRMRD::Dataset dataset("csm_testdata.h5", "dataset");
		//ISMRMRD::Dataset dataset("csm.h5", get_date_time_string().c_str());
		dataset.appendImage("csm", img);
		mtx.unlock();
		//iw.write(dataset);
		delete ptr_img;
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}
#endif

int test10()
{
	Mutex mtx;
	mtx.lock();
	//ISMRMRD::Dataset input("opismrmrd_csm.h5", "coil_maps");
	//ISMRMRD::Dataset output("csm_opismrmrd.h5", "dataset");
	ISMRMRD::Dataset input("opismrmrd_with_csm.h5", "coil_maps");
	ISMRMRD::Dataset output("csm1.h5", "dataset");
	mtx.unlock();

	int ni = input.getNumberOfImages("image_1000");
	std::cout << ni << std::endl;
	int mi = input.getNumberOfImages("image_2000");
	std::cout << mi << std::endl;

	ISMRMRD::Image<float> re;
	ISMRMRD::Image<float> im;
	CFImage img;

	for (int i = 0; i < ni; i++) {
		mtx.lock();
		input.readImage("image_1000", i, re);
		input.readImage("image_2000", i, im);
		mtx.unlock();
		int nx = re.getMatrixSizeX();
		int ny = re.getMatrixSizeY();
		int nz = re.getMatrixSizeZ();
		int nc = re.getNumberOfChannels();
		if (i == 0) {
			std::cout << nx << std::endl;
			std::cout << ny << std::endl;
			std::cout << nz << std::endl;
			std::cout << nc << std::endl;
		}
		ISMRMRD::Image<complex_float_t>* ptr_img =
			new ISMRMRD::Image<complex_float_t>(nx, ny, nz, nc);
		//ImageWrap iw(ISMRMRD::ISMRMRD_CXFLOAT, ptr_img);
		ISMRMRD::Image<complex_float_t>& img = *ptr_img;
		for (int ic = 0; ic < nc; ic++)
			for (int iz = 0; iz < nz; iz++)
				for (int iy = 0; iy < ny; iy++)
					for (int ix = 0; ix < nx; ix++)
						img(ix, iy, iz, ic) = 
							complex_float_t(re(ix, iy, iz, ic), im(ix, iy, iz, ic));
		mtx.lock();
		output.appendImage("csm", img);
		mtx.unlock();
		delete ptr_img;
	}

	//for (int i = 0; i < ni; i++) {
	//	input.readImage("image_1000", i, re);
	//	int nx = re.getMatrixSizeX();
	//	int ny = re.getMatrixSizeY();
	//	int nz = re.getMatrixSizeZ();
	//	int nc = re.getNumberOfChannels();
	//	int ns = re.getSlice();
	//	//std::cout << nx << std::endl;
	//	//std::cout << ny << std::endl;
	//	//std::cout << nz << std::endl;
	//	//std::cout << nc << std::endl;
	//	std::cout << ' ' << ns;
	//}
	//std::cout << std::endl;
	//for (int i = 0; i < ni; i++) {
	//	input.readImage("image_2000", i, im);
	//	int nx = im.getMatrixSizeX();
	//	int ny = im.getMatrixSizeY();
	//	int nz = im.getMatrixSizeZ();
	//	int nc = im.getNumberOfChannels();
	//	int ns = im.getSlice();
	//	//std::cout << nx << std::endl;
	//	//std::cout << ny << std::endl;
	//	//std::cout << nz << std::endl;
	//	//std::cout << nc << std::endl;
	//	std::cout << ' ' << ns;
	//}
	//std::cout << std::endl;

	return 0;
}

int test11(const char* in_file)
{
	//CoilSensitivitiesAsImages csms(in_file);
	//int n = csms.items();
	//void* ptr_csms = cGT_CoilSensitivitiesFromFile(in_file);
	void* ptr_csms = cGT_CoilSensitivities(in_file);
	void* ptr_n = cGT_dataItems(ptr_csms);
	int n = intDataFromHandle(ptr_n);

	int dim[4];
	//csms.get_dim(0, dim);
	cGT_getCSMDimensions(ptr_csms, 0, (size_t)dim);
	std::cout << n << std::endl;
	int nx = dim[0];
	int ny = dim[1];
	int nz = dim[2];
	int nc = dim[3];

	//const CFImage& csm = csms(0);
	//int nx = csm.getMatrixSizeX();
	//int ny = csm.getMatrixSizeY();
	//int nz = csm.getMatrixSizeZ();
	//int nc = csm.getNumberOfChannels();

	std::cout << nx << std::endl;
	std::cout << ny << std::endl;
	std::cout << nz << std::endl;
	std::cout << nc << std::endl;

	long long int size = nx;
	size *= ny;
	size *= nz;
	size *= nc;
	double* re = new double[size];
	double* im = new double[size];
	double* v = new double[size];

	//csms.get_data(1, re, im);
	cGT_getCSMData(ptr_csms, 1, (size_t)re, (size_t)im);
	cGT_getCSMDataAbs(ptr_csms, 1, (size_t)v);

	std::cout << re[0] << std::endl;
	std::cout << im[0] << std::endl;
	std::cout << v[0] << std::endl;

	std::cout << re[1] << std::endl;
	std::cout << im[1] << std::endl;
	std::cout << v[1] << std::endl;

	delete[] re;
	delete[] im;
	delete[] v;

	deleteObject(ptr_csms);
	deleteObject(ptr_n);

	return 0;
}

#include <algorithm>
#include <functional>
#include <array>
#include <iostream>

class pair {
public:
	pair(int f, int s) : first(f), second(s) {}
	int first;
	int second;
};

int test12()
{
	std::array<int, 10> s = { 5, 7, 4, 2, 8, 6, 1, 9, 0, 3 };

	//// sort using the default operator<
	//std::sort(s.begin(), s.end());
	//for (int a : s) {
	//	std::cout << a << " ";
	//}
	//std::cout << '\n';

	std::vector<pair> p;
	//= {
	//	pair(0,5), pair(1,7), pair(2,4), pair(3,2), pair(4,8), pair(5,6), 
	//	pair(6,1), pair(7,9), pair(8,0), pair(9,3) };

	for (int i = 0; i < 10; i++)
		p.push_back(pair(i, s[i]));

	std::sort(p.begin(), p.end(), [](pair a, pair b) {return b.second > a.second;});

	for (pair a : p) {
		std::cout << a.first << ' ';
	}
	std::cout << '\n';

	return 0;
}

int test13(const char* in_file)
{
	std::cout << "Gadgetron ISMRMRD client" << std::endl;
	std::cout << "  -- hdf5 file  in   :      " << in_file << std::endl;

	try {
		boost::shared_ptr<AcquisitionsContainer> sptr_input(new AcquisitionsFile(in_file));
		AcquisitionsContainer& input = *sptr_input;
		input.order();
		for (int i = 0; i < input.number(); i++)
			std::cout << input.index(i) << ' ';
		std::cout << '\n';
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}
