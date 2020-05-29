/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 Rutherford Appleton Laboratory STFC

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
\ingroup Gadgetron Extensions
\brief Auxiliary functions for MR related C++ tests.

\author Johannes Mayer
\author CCP PETMR
*/

#pragma once

#include "sirf/Gadgetron/gadgetron_x.h"
#include "sirf/Gadgetron/gadget_lib.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_image_wrap.h"

namespace sirf{


void preprocess_acquisition_data(MRAcquisitionData& ad);

void write_cfimage_to_raw(const std::string& fname_prefix, const CFImage& img);
void write_cfimage_to_raw(const std::string& fname_prefix, const ImageWrap& iw);

} // END NAMESPACE
