/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 Rutherford Appleton Laboratory STFC
Copyright 2020 University College London

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

/*!
\file
\ingroup SIRF C++ Synergistic Utilities
\brief C++ synergistic utilities.

\author Richard Brown
\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/STIR/stir_data_containers.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Syn/syn_utilities.h"

using namespace sirf;

ImageDataWrap::ImageDataWrap(const std::string &filename, const std::string &engine, bool verbose)
{
    if (strcmp(engine.c_str(), "Reg") == 0) {
        std::shared_ptr<NiftiImageData<float> > nifti_sptr = 
		std::make_shared<NiftiImageData3D<float> >(filename);
        if (verbose) nifti_sptr->print_header();
        img_sptr_ = nifti_sptr;
    }
    else if (strcmp(engine.c_str(), "STIR") == 0) {
        img_sptr_ = std::make_shared<STIRImageData>(filename);
    }
    else if (strcmp(engine.c_str(), "Gadgetron") == 0) {
        std::shared_ptr<GadgetronImagesVector> gadgetron_sptr(new GadgetronImagesVector);
		gadgetron_sptr->read(filename);
        if (verbose) gadgetron_sptr->print_header(0);
        img_sptr_ = gadgetron_sptr;
    }
    else
        throw std::runtime_error("unknown engine - " + engine + ".\n");

    // If verbose print geom info
	if (verbose) {
		std::shared_ptr<const VoxelisedGeometricalInfo3D > gi_sptr =
			img_sptr_->get_geom_info_sptr();
		if (gi_sptr.get())
			gi_sptr->print_info();
	}
}