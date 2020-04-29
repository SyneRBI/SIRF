// author: Johannes Mayer
// date: 2020.02.27


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
