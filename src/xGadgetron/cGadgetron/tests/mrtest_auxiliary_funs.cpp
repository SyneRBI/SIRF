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

#include "mrtest_auxiliary_funs.h"


#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

#include <ismrmrd/ismrmrd.h>

void sirf::preprocess_acquisition_data(MRAcquisitionData& ad)
{
    std::cout << "Pre-processing Acquisition Data" << std::endl;

    sirf::AcquisitionsProcessor preprocessing_chain;

    auto sptr_noise_gadget = std::make_shared<Gadget>(NoiseAdjustGadget());
    auto sptr_ro_overs_gadget = std::make_shared<Gadget>(RemoveROOversamplingGadget());
    auto sptr_asymmecho_gadget = std::make_shared<Gadget>(AsymmetricEchoAdjustROGadget());

    preprocessing_chain.add_gadget("dummy1", sptr_noise_gadget);
    preprocessing_chain.add_gadget("dummy2", sptr_asymmecho_gadget);
    preprocessing_chain.add_gadget("dummy3", sptr_ro_overs_gadget);

    preprocessing_chain.process(ad);
    auto sptr_preproc_ad =preprocessing_chain.get_output();

    ISMRMRD::Acquisition acq;
    for(int i=0; i<sptr_preproc_ad->number(); ++i)
    {
        sptr_preproc_ad->get_acquisition(i, acq);
        ad.set_acquisition(i, acq);
    }
    ad.set_acquisitions_info( sptr_preproc_ad->acquisitions_info());

}

void sirf::write_cfimage_to_raw(const std::string& fname_prefix, const CFImage& img)
{


    std::stringstream fname_out;
    fname_out << fname_prefix;
    fname_out << "_" << img.getMatrixSizeX();
    fname_out << "x" << img.getMatrixSizeY();
    fname_out << "x" << img.getMatrixSizeZ() * img.getNumberOfChannels();
    fname_out << ".raw";

    std::cout << "Writing " << fname_out.str() << std::endl;

    std::vector<float> dat;

    for(size_t i=0; i<img.getNumberOfDataElements(); ++i)
        dat.push_back( std::abs( *(img.getDataPtr() + i )));


    std::ofstream myfile( fname_out.str(), std::ios::out | std::ios::binary);
    myfile.write((char*)(&dat[0]), sizeof(float)*dat.size());

}


void sirf::write_cfimage_to_raw(const std::string& fname_prefix, const ImageWrap& iw)
{
    const void* vptr_img = iw.ptr_image();
    const CFImage* ptr_img = static_cast<const CFImage*>(vptr_img);

    sirf::write_cfimage_to_raw(fname_prefix, *ptr_img);
}
