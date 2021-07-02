/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2018 Rutherford Appleton Laboratory STFC

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

/*!
\file
\ingroup PET

SIRF C interface test.

\author Evgueni Ovtchinnikov
\author SyneRBI
*/
#include <string>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"

#include "sirf/common/csirf.h"
#include "sirf/STIR/cstir.h"
#include "sirf/STIR/stir_types.h"

#include "getenv.h"
#include "handle.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

void* TMP_HANDLE;

#define GET_FLOAT(V, F) \
	TMP_HANDLE = F; \
	if (execution_status(TMP_HANDLE)) break; \
	V = floatDataFromHandle(TMP_HANDLE); \
	deleteDataHandle(TMP_HANDLE)

int test2()
{
	std::cout << "running test2.cpp...\n";
	std::string filename;
	int dim[10];
	float at_value = 0.05f*0;
	float bt_value = 0.1f*0;
	float s, t;
	void* handle = 0;
	void* matrix = 0;
	void* image = 0;
	void* img = 0;
	void* am = 0;
	void* ad = 0;
	void* fd = 0;
	void* at = 0;
	void* bt = 0;
	void* nd = 0;
	void* norm = 0;
	void* prior = 0;
	void* obj_fun = 0;
	void* filter = 0;
	void* recon = 0;
	void* diff = 0;
	void* sm = 0;
	void* sma = 0;
	void* smn = 0;
	void* ai = 0;

	int status = 1;
	try {
		std::string SIRF_path = sirf::getenv("SIRF_PATH");
		if (SIRF_path.length() < 1) {
			std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
			return 1;
		}
		std::string path = SIRF_path + "/data/examples/PET/";

		TextWriter w;
		openChannel(0, &w); // suppress STIR info output

		for (;;) {

			HANDLE(matrix, cSTIR_newObject("RayTracingMatrix"));
			CALL(cSTIR_setParameter
			(matrix, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2)));

			filename = path + "mMR/mMR_template_span11_small.hs";
			//filename = path + "my_forward_projection.hs";
			//std::cout << "reading data from " << filename << "...";
			//BUG: fails if storage scheme is "memory"!
			HANDLE(ad, cSTIR_objectFromFile("AcquisitionData", filename.c_str()));
			//std::cout << "ok\n";
			cSTIR_getAcquisitionDataDimensions(ad, (size_t)&dim[0]);
			std::cout << "acquisition data dimensions: "
				<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

			//cSTIR_setAcquisitionDataStorageScheme("memory");

			HANDLE(image, cSTIR_imageFromAcquisitionData(ad));
			cSTIR_getImageDimensions(image, (size_t)&dim[0]);
			std::cout << "image dimensions: "
				<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

			HANDLE(at, cSTIR_acquisitionDataFromTemplate(ad));
			cSTIR_getAcquisitionDataDimensions(at, (size_t)&dim[0]);
			cSTIR_fillAcquisitionData(at, at_value);
			HANDLE(bt, cSTIR_acquisitionDataFromTemplate(ad));
			cSTIR_getAcquisitionDataDimensions(bt, (size_t)&dim[0]);
			cSTIR_fillAcquisitionData(bt, bt_value);
			HANDLE(nd, cSTIR_acquisitionDataFromTemplate(ad));
			cSTIR_getAcquisitionDataDimensions(nd, (size_t)&dim[0]);
			cSTIR_fillAcquisitionData(nd, 2.0);

			HANDLE(am, cSTIR_newObject("AcqModUsingMatrix"));
			CALL(cSTIR_setParameter(am, "AcquisitionModel", "additive_term", at));
			CALL(cSTIR_setParameter(am, "AcqModUsingMatrix", "matrix", matrix));
			CALL(cSTIR_setupAcquisitionModel(am, ad, image));

			filename = path + "mMR/mu_map.hv";
			HANDLE(ai, cSTIR_objectFromFile("Image", filename.c_str()));
			HANDLE(sma, cSTIR_createPETAttenuationModel(ai, am));
			CALL(cSTIR_setupAcquisitionSensitivityModel(sma, ad));
			HANDLE(smn, cSTIR_createPETAcquisitionSensitivityModel(nd, "s"));
			CALL(cSTIR_fillAcquisitionData(nd, 1.0));
			CALL(cSTIR_applyAcquisitionSensitivityModel(sma, nd, "unnormalise"));
			deleteDataHandle(sma);
			HANDLE(sma, cSTIR_createPETAcquisitionSensitivityModel(nd, "s"));
			HANDLE(sm, cSTIR_chainPETAcquisitionSensitivityModels(smn, sma));
			CALL(cSTIR_setupAcquisitionSensitivityModel(sm, ad));
			CALL(cSTIR_setParameter(am, "AcquisitionModel", "asm", sm));

			int num_subsets = 9; // 8;
			std::cout << "projecting...\n";
			HANDLE(fd, cSTIR_acquisitionModelFwd(am, image, 0, num_subsets));
			cSTIR_getAcquisitionDataDimensions(fd, (size_t)&dim[0]);
			std::cout << "simulated acquisition data dimensions: "
				<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

			std::cout << "backprojecting...\n";
			HANDLE(img, cSTIR_acquisitionModelBwd(am, fd, 0, num_subsets));
			cSTIR_getImageDimensions(img, (size_t)&dim[0]);
			std::cout << "backprojected image dimensions: "
				<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

			HANDLE(prior, cSTIR_newObject("QuadraticPrior"));

			std::string obj_fun_name
			("PoissonLogLikelihoodWithLinearModelForMeanAndProjData");
			HANDLE(obj_fun, cSTIR_newObject(obj_fun_name.c_str()));
			CALL(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "acquisition_model", am));
			handle = charDataHandle("true");
			CALL(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "zero_seg0_end_planes", handle));
			deleteDataHandle(handle);
			CALL(cSTIR_setParameter
			(obj_fun, "GeneralisedObjectiveFunction", "prior", prior));

			HANDLE(filter, cSTIR_newObject("TruncateToCylindricalFOVImageProcessor"));

			int num_subiterations = 2;
			HANDLE(recon, cSTIR_objectFromFile("OSMAPOSLReconstruction", ""));
			handle = charDataHandle("reconstructedImage");
			CALL(cSTIR_setParameter
			(recon, "Reconstruction", "output_filename_prefix", handle));
			deleteDataHandle(handle);
			handle = intDataHandle(num_subsets);
			CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "num_subsets", handle));
			deleteDataHandle(handle);
			handle = intDataHandle(num_subiterations);
			CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "num_subiterations", handle));
			CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "save_interval", handle));
			deleteDataHandle(handle);
			handle = intDataHandle(1);
			CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "inter_iteration_filter_interval",
				handle));
			deleteDataHandle(handle);
			CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "objective_function", obj_fun));
			CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "inter_iteration_filter_type", filter));
			handle = charDataHandle("multiplicative");
			CALL(cSTIR_setParameter(recon, "OSMAPOSL", "MAP_model", handle));
			deleteDataHandle(handle);
			CALL(cSTIR_setParameter(recon, "Reconstruction", "input_data", fd));
			std::cout << "setting up the reconstructor, please wait...";
			CALL(cSTIR_setupReconstruction(recon, image));
			std::cout << "ok\n";

			for (int iter = 0; iter < num_subiterations; iter++) {
				std::cout << "iteration " << iter << '\n';
				cSTIR_updateReconstruction(recon, image);
			}

			status = 0;
			break;
		}

		deleteDataHandle(image);
		deleteDataHandle(img);
		deleteDataHandle(ad);
		deleteDataHandle(am);
		deleteDataHandle(sm);
		deleteDataHandle(sma);
		deleteDataHandle(smn);
		deleteDataHandle(matrix);
		deleteDataHandle(at);
		deleteDataHandle(bt);
		deleteDataHandle(nd);
		deleteDataHandle(fd);
		deleteDataHandle(recon);
		deleteDataHandle(filter);
		deleteDataHandle(prior);
		deleteDataHandle(obj_fun);
		closeChannel(0, &w);
		std::cout << "done with test2.cpp...\n";
	}
	catch (const std::exception &error) {
		std::cerr << "\nException thrown:\n\t" << error.what() << "\n\n";
		return 1;
	}
	return status;
}

