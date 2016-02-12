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

int test1(
const char* host_name,
const char* port,
unsigned int timeout_ms,
const char* in_filename,
const char* hdf5_in_group,
const char* out_filename,
const char* hdf5_out_group,
const char* config_file);

int test2(
const char* host_name,
const char* port,
unsigned int timeout_ms,
const char* in_filename,
const char* hdf5_in_group,
const char* out_filename,
const char* hdf5_out_group,
const char* config_file);

int test3(
const char* host_name,
const char* port,
unsigned int timeout_ms,
const char* in_filename,
const char* hdf5_in_group,
const char* out_filename,
const char* hdf5_out_group);
//const char* config_file);

int test4();

int test5(
	const char* host,
	const char* port,
	unsigned int timeout,
	const char* in_file,
	const char* in_group,
	const char* out_file,
	const char* out_group);

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

namespace po = boost::program_options;
using boost::asio::ip::tcp;

int main(int argc, char **argv)
{
	//return test4();

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

	return test7(
		host_name.c_str(),
		port.c_str(),
		timeout_ms,
		in_filename.c_str(),
		hdf5_in_group.c_str(),
		"out7.h5",
		hdf5_out_group.c_str());

	return test6(
		host_name.c_str(),
		port.c_str(),
		timeout_ms,
		in_filename.c_str(),
		hdf5_in_group.c_str(),
		"out6.h5",
		hdf5_out_group.c_str());

	return test5(
		host_name.c_str(),
		port.c_str(),
		timeout_ms,
		in_filename.c_str(),
		hdf5_in_group.c_str(),
		"out5.h5",
		hdf5_out_group.c_str());

	return test3(
		host_name.c_str(),
		port.c_str(),
		timeout_ms,
		in_filename.c_str(),
		hdf5_in_group.c_str(),
		"out3.h5",
		hdf5_out_group.c_str());

	if (vm.count("config-local")) {
		std::ifstream t(config_file_local.c_str());
		if (t) {
			//Read in the file.
			config_xml_local = std::string((std::istreambuf_iterator<char>(t)),
				std::istreambuf_iterator<char>());
			return test3(
				host_name.c_str(), 
				port.c_str(),
				 timeout_ms,
				in_filename.c_str(),
				hdf5_in_group.c_str(),
				"out3.h5",
				 hdf5_out_group.c_str());
				 //config_xml_local.c_str());
		}
		else {
			std::cout << "Unable to read local xml configuration: " << config_file_local << std::endl;
			return -1;
		}
	}

	//TODO:
	// Add check to see if input file exists

	return test2(
		host_name.c_str(), 
		port.c_str(),
		 timeout_ms,
		in_filename.c_str(),
		hdf5_in_group.c_str(),
		"out2.h5",
		 hdf5_out_group.c_str(),
		 config_file.c_str());

	return 0;
}

int test1(
const char* host,
const char* port,
unsigned int timeout,
const char* in_file,
const char* in_group,
const char* out_file,
const char* out_group,
const char* config_file)
{
	std::cout << "Gadgetron ISMRMRD client" << std::endl;
	std::cout << "  -- host            :      " << host << std::endl;
	std::cout << "  -- port            :      " << port << std::endl;
	std::cout << "  -- hdf5 file  in   :      " << in_file << std::endl;
	std::cout << "  -- hdf5 group in   :      " << in_group << std::endl;
	std::cout << "  -- conf            :      " << config_file << std::endl;
	std::cout << "  -- hdf5 file out   :      " << out_file << std::endl;
	std::cout << "  -- hdf5 group out  :      " << out_group << std::endl;
	
	try {
		void* h_data = cGT_ISMRMRDatasetFromFile(in_file, in_group);
		void* h_head = cGT_newObject("string");
		void* h_conn = cGT_newObject("GTConnector");

		cGT_readISMRMRDatasetHeader(h_data, h_head);
		cGT_setConnectionTimeout(h_conn, timeout);
		cGT_registerHDFReceiver(h_conn, out_file, out_group);
		cGT_connect(h_conn, host, port);
		cGT_sendConfigFile(h_conn, config_file);	
		cGT_sendParameters(h_conn, h_head);
		cGT_sendAcquisitions(h_conn, h_data);
		cGT_disconnect(h_conn);

		deleteObject(h_data);
		deleteObject(h_head);
		deleteObject(h_conn);
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}

int test2(
const char* host,
const char* port,
unsigned int timeout,
const char* in_file,
const char* in_group,
const char* out_file,
const char* out_group,
const char* config_file)
{
	std::cout << "Gadgetron ISMRMRD client" << std::endl;
	std::cout << "  -- host            :      " << host << std::endl;
	std::cout << "  -- port            :      " << port << std::endl;
	std::cout << "  -- hdf5 file  in   :      " << in_file << std::endl;
	std::cout << "  -- hdf5 group in   :      " << in_group << std::endl;
	std::cout << "  -- conf            :      " << config_file << std::endl;
	std::cout << "  -- hdf5 file out   :      " << out_file << std::endl;
	std::cout << "  -- hdf5 group out  :      " << out_group << std::endl;
	
	try {
		void* h_data = cGT_ISMRMRDatasetFromFile(in_file, in_group);
		void* h_head = cGT_newObject("string");
		void* h_conn = cGT_newObject("GTConnector");
		void* h_imgs = cGT_newObject("ImagesList");

		cGT_readISMRMRDatasetHeader(h_data, h_head);
		cGT_setConnectionTimeout(h_conn, timeout);
		cGT_registerImagesReceiver(h_conn, h_imgs);
		cGT_connect(h_conn, host, port);
		cGT_sendConfigFile(h_conn, config_file);	
		cGT_sendParameters(h_conn, h_head);
		cGT_sendAcquisitions(h_conn, h_data);
		cGT_disconnect(h_conn);

		cGT_writeImages(h_imgs, out_file, out_group);
		//cGT_writeImages(h_imgs, h_conn, out_file, out_group);

		int im_num;
		std::cout << "image: ";
		std::cin >> im_num;
		int dim[3];
		size_t ptr_dim = (size_t)dim;
		cGT_getImageDimensions(h_imgs, im_num, ptr_dim);
		std::cout << "image sizes: " << dim[0] << 'x' << dim[1] << 'x' << dim[2] << std::endl;
		size_t n = dim[0];
		n *= dim[1];
		n *= dim[2];
		double* ptr = new double[n];
		size_t ptr_data = (size_t)ptr;
		cGT_getImageDataAsDoubleArray(h_imgs, im_num, ptr_data);
		int sx = dim[0];
		int sy = dim[1];
		int sz = dim[2];
		std::cout << ptr[0] << '\t' << ptr[1] << std::endl;
		std::cout << ptr[sx] << '\t' << ptr[sx + 1] << std::endl;
		std::cout << ptr[2*sx] << '\t' << ptr[2*sx + 1] << std::endl;
		delete[] ptr;

		deleteObject(h_data);
		deleteObject(h_head);
		deleteObject(h_conn);
		deleteObject(h_imgs);
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}

int test3(
const char* host,
const char* port,
unsigned int timeout,
const char* in_file,
const char* in_group,
const char* out_file,
const char* out_group)
//const char* config_file)
{
	std::cout << "Gadgetron ISMRMRD client" << std::endl;
	std::cout << "  -- host            :      " << host << std::endl;
	std::cout << "  -- port            :      " << port << std::endl;
	std::cout << "  -- hdf5 file  in   :      " << in_file << std::endl;
	std::cout << "  -- hdf5 group in   :      " << in_group << std::endl;
	//std::cout << "  -- conf            :      " << config_file << std::endl;
	std::cout << "  -- hdf5 file out   :      " << out_file << std::endl;
	std::cout << "  -- hdf5 group out  :      " << out_group << std::endl;
	
	try {
		void* h_data = cGT_ISMRMRDatasetFromFile(in_file, in_group);
		void* h_head = cGT_newObject("string");
		void* h_conn = cGT_newObject("GTConnector");
		void* h_imgs = cGT_newObject("ImagesList");

		cGT_readISMRMRDatasetHeader(h_data, h_head);
		cGT_setConnectionTimeout(h_conn, timeout);
		cGT_registerImagesReceiver(h_conn, h_imgs);
		//cGT_registerHDFReceiver(h_conn, out_file, out_group);
		cGT_connect(h_conn, host, port);
		void* h_gc = cGT_newObject("GadgetChain");
		void* h_ar = cGT_newObject("GadgetIsmrmrdAcquisitionMessageReader");
		void* h_iw = cGT_newObject("MRIImageWriter");
		void* h_ro = cGT_newObject("RemoveROOversamplingGadget");
		void* h_aat = cGT_newObject("AcquisitionAccumulateTriggerGadget");
		void* h_bb = cGT_newObject("BucketToBufferGadget");
		void* h_sr = cGT_newObject("SimpleReconGadget");
		void* h_ias = cGT_newObject("ImageArraySplitGadget");
		void* h_e = cGT_newObject("ExtractGadget");
		void* h_if = cGT_newObject("ImageFinishGadget");

		void* status; 
		status = cGT_setGadgetProperty(h_aat, "trigger_dimension", "repetition");
		if (executionStatus(status)) {
			std::cout << "exception thrown" << std::endl;
			exit(1);
		}
		status = cGT_setGadgetProperty(h_bb, "n_dimension", "");
		if (executionStatus(status)) {
			std::cout << "exception thrown" << std::endl;
			exit(1);
		}

		cGT_addReader(h_gc, "reader1", h_ar);
		cGT_addWriter(h_gc, "writer1", h_iw);
		cGT_addGadget(h_gc, "gadget1", h_ro);
		cGT_addGadget(h_gc, "gadget2", h_aat);
		cGT_addGadget(h_gc, "gadget3", h_bb);
		cGT_addGadget(h_gc, "gadget4", h_sr);
		cGT_addGadget(h_gc, "gadget5", h_ias);
		//cGT_addGadget(h_gc, "gadget6", h_e);
		//cGT_addGadget(h_gc, "endgadget", h_if);
		cGT_setEndGadget(h_gc, h_if);

		cGT_configGadgetChain(h_conn, h_gc);
		
		cGT_sendParameters(h_conn, h_head);
		cGT_sendAcquisitions(h_conn, h_data);
		cGT_disconnect(h_conn);
		
		void* h_gc2 = cGT_newObject("GadgetChain");
		void* h_imgs2 = cGT_newObject("ImagesList");
		void* h_ir = cGT_newObject("MRIImageReader");
		cGT_registerImagesReceiver(h_conn, h_imgs2);
		cGT_addReader(h_gc2, "reader2", h_ir);
		cGT_addWriter(h_gc2, "writer1", h_iw);
		//cGT_addGadget(h_gc, "gadget5", h_ias);
		cGT_addGadget(h_gc2, "gadget6", h_e);
		//cGT_addGadget(h_gc2, "endgadget", h_if);
		cGT_setEndGadget(h_gc2, h_if);

		cGT_connect(h_conn, host, port);
		cGT_configGadgetChain(h_conn, h_gc2);
		cGT_sendImages(h_conn, h_imgs);
		cGT_disconnect(h_conn);

		cGT_writeImages(h_imgs2, out_file, out_group);
		//cGT_writeImages(h_imgs2, h_conn, out_file, out_group);

		//cGT_writeImages(h_imgs, h_conn, out_file, out_group);

		deleteObject(h_data);
		deleteObject(h_head);
		deleteObject(h_conn);
		deleteObject(h_imgs);
		deleteObject(h_gc);
		deleteObject(h_ar);
		deleteObject(h_iw);
		deleteObject(h_ro);
		deleteObject(h_aat);
		deleteObject(h_bb);
		deleteObject(h_sr);
		deleteObject(h_ias);
		deleteObject(h_e);
		deleteObject(h_if);

		deleteObject(h_imgs2);
		deleteObject(h_gc2);
		deleteObject(h_ir);
		//deleteObject(h_ir);
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}
	catch (...) {
		std::cout << "unhandled exception thrown" << std::endl;
	}

	return 0;
}

int test4() {

/*
	GadgetChain gc;
	NEW(IsmrmrdAcqMsgReader, ptr_r);
	gc.add_reader("reader1", ptr_r );
	gc.add_writer("writer1", new IsmrmrdImgMsgWriter);
	gc.add_gadget("gadget1", new RemoveOversamplingGadget);
	gc.add_gadget("gadget2", new AcqAccTrigGadget);
	gc.add_gadget("gadget3", new BucketToBuffGadget);
	gc.add_gadget("gadget4", new SimpleReconstructionGadget);
	gc.add_gadget("gadget5", new ImgArrSplitGadget);
	gc.add_gadget("gadget6", new ExtGadget);
	gc.add_gadget("gadget7", new ImgFinishGadget);
	
	std::cout << gc.xml() << std::endl;
	
	delete ptr_r;
*/

	try {
		void* h_conn = cGT_newObject("GTConnector");
		void* h_gc = cGT_newObject("GadgetChain");
		void* h_r = cGT_newObject("IsmrmrdAcqMsgReader");
		void* h_w = cGT_newObject("IsmrmrdImgMsgWriter");
		void* h_ro = cGT_newObject("RemoveOversamplingGadget");
		void* h_aat = cGT_newObject("AcqAccTrigGadget");
		void* h_bb = cGT_newObject("BucketToBuffGadget");
		void* h_sr = cGT_newObject("SimpleReconstructionGadget");
		void* h_ias = cGT_newObject("ImgArrSplitGadget");
		void* h_e = cGT_newObject("ExtGadget");
		void* h_if = cGT_newObject("ImgFinishGadget");
		
		cGT_addReader(h_gc, "reader1", h_r);
		cGT_addWriter(h_gc, "writer1", h_w);
		cGT_addGadget(h_gc, "gadget1", h_ro);
		cGT_addGadget(h_gc, "gadget2", h_aat);
		cGT_addGadget(h_gc, "gadget3", h_bb);
		cGT_addGadget(h_gc, "gadget4", h_sr);
		cGT_addGadget(h_gc, "gadget5", h_ias);
		cGT_addGadget(h_gc, "gadget6", h_e);
		cGT_addGadget(h_gc, "gadget7", h_if);

		cGT_configGadgetChain(h_conn, h_gc);
		
		deleteObject(h_gc);
		deleteObject(h_r);
		deleteObject(h_w);
		deleteObject(h_ro);
		deleteObject(h_aat);
		deleteObject(h_bb);
		deleteObject(h_sr);
		deleteObject(h_ias);
		deleteObject(h_e);
		deleteObject(h_if);
		deleteObject(h_conn);
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}
	catch (...) {
		std::cout << "unhandled exception thrown" << std::endl;
		return -1;
	}

	return 0;
}

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
		ISMRMRD::Dataset input(in_file, in_group, false);
		std::cout << "ok" << std::endl;

		SPTR(aGadget, ro, RemoveOversamplingGadget);
		SPTR(aGadget, aat, AcqAccTrigGadget);
		SPTR(aGadget, bb, BucketToBuffGadget);
		SPTR(aGadget, sr, SimpleReconstructionGadget);
		SPTR(aGadget, ias, ImgArrSplitGadget);
		SPTR(aGadget, ext, ExtGadget);
		SPTR(aGadget, fin, ImgFinishGadget);
		std::cout << "ok" << std::endl;

		{
			MRIReconstruction recon;
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

			ImagesProcessor proc;
			proc.add_gadget("g1", ext);
			//proc.add_gadget("g2", fin);
			proc.process(imgs);
			ImagesList& images = (ImagesList&)*proc.get_output();

			std::cout << images.number() << std::endl;

			//images.write(out_file, out_group, GTConnector());
			images.write(out_file, out_group);
			//std::cout << "ok" << std::endl;

		}
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}

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

		void* h_recon = cGT_newObject("MRIReconstruction");
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

		std::cout << cGT_numImages(h_imgs) << std::endl;

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

	try {
		void* h_input = cGT_ISMRMRDAcquisitionsFromFile(in_file);

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
		void* h_proc = cGT_acquisitionsProcessor(acq_file.c_str());
		cGT_addGadget(h_proc, "g1", h_ro);

		void* h_output = cGT_processAcquisitions(h_proc, h_input);
		if (executionStatus(h_output)) {
			std::cout << "exception thrown" << std::endl;
			exit(1);
		}

		void* h_recon = cGT_newObject("MRIReconstruction");
		cGT_addGadget(h_recon, "g1", h_sr);
		cGT_addGadget(h_recon, "g2", h_e);

		void* h_images = cGT_reconstructImages(h_recon, h_output);

		cGT_writeImages(h_images, out_file, out_group);

		deleteObject(h_ro);
		deleteObject(h_sr);
		deleteObject(h_e);
		deleteObject(h_recon);
		deleteObject(h_input);
		deleteObject(h_output);
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}

