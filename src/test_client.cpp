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

namespace po = boost::program_options;
using boost::asio::ip::tcp;

int main(int argc, char **argv)
{
	return test4();

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
		void* h_head = newObject("string");
		void* h_conn = newObject("GTConnector");

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
		void* h_head = newObject("string");
		void* h_conn = newObject("GTConnector");
		void* h_imgs = newObject("ImagesList");

		cGT_readISMRMRDatasetHeader(h_data, h_head);
		cGT_setConnectionTimeout(h_conn, timeout);
		cGT_registerImagesReceiver(h_conn, h_imgs);
		cGT_connect(h_conn, host, port);
		cGT_sendConfigFile(h_conn, config_file);	
		cGT_sendParameters(h_conn, h_head);
		cGT_sendAcquisitions(h_conn, h_data);
		cGT_disconnect(h_conn);

		cGT_writeImages(h_imgs, h_conn, out_file, out_group);

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
		void* h_head = newObject("string");
		void* h_conn = newObject("GTConnector");

		cGT_readISMRMRDatasetHeader(h_data, h_head);
		cGT_setConnectionTimeout(h_conn, timeout);
		cGT_registerHDFReceiver(h_conn, out_file, out_group);
		cGT_connect(h_conn, host, port);
/*
		GadgetChain gc;
		gc.add_reader("reader1", new IsmrmrdAcqMsgReader);
		gc.add_writer("writer1", new IsmrmrdImgMsgWriter);
		gc.add_gadget("gadget1", new RemoveOversamplingGadget);
		gc.add_gadget("gadget2", new AcqAccTrigGadget);
		gc.add_gadget("gadget3", new BucketToBuffGadget);
		gc.add_gadget("gadget4", new SimpleReconstructionGadget);
		gc.add_gadget("gadget5", new ImgArrSplitGadget);
		gc.add_gadget("gadget6", new ExtGadget);
		gc.add_gadget("gadget7", new ImgFinishGadget);
		std::string config_file = gc.xml();
		cGT_sendConfigScript(h_conn, config_file.c_str());	
*/
		void* h_gc = newObject("GadgetChain");
		void* h_r = newObject("IsmrmrdAcqMsgReader");
		void* h_w = newObject("IsmrmrdImgMsgWriter");
		void* h_ro = newObject("RemoveOversamplingGadget");
		void* h_aat = newObject("AcqAccTrigGadget");
		void* h_bb = newObject("BucketToBuffGadget");
		void* h_sr = newObject("SimpleReconstructionGadget");
		void* h_ias = newObject("ImgArrSplitGadget");
		void* h_e = newObject("ExtGadget");
		void* h_if = newObject("ImgFinishGadget");
		
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
		
		cGT_sendParameters(h_conn, h_head);
		cGT_sendAcquisitions(h_conn, h_data);
		cGT_disconnect(h_conn);

		deleteObject(h_data);
		deleteObject(h_head);
		deleteObject(h_conn);
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
	}
	catch (std::exception& ex) {
		std::cout << "Error caught: " << ex.what() << std::endl;
		return -1;
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
		void* h_conn = newObject("GTConnector");
		void* h_gc = newObject("GadgetChain");
		void* h_r = newObject("IsmrmrdAcqMsgReader");
		void* h_w = newObject("IsmrmrdImgMsgWriter");
		void* h_ro = newObject("RemoveOversamplingGadget");
		void* h_aat = newObject("AcqAccTrigGadget");
		void* h_bb = newObject("BucketToBuffGadget");
		void* h_sr = newObject("SimpleReconstructionGadget");
		void* h_ias = newObject("ImgArrSplitGadget");
		void* h_e = newObject("ExtGadget");
		void* h_if = newObject("ImgFinishGadget");
		
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
