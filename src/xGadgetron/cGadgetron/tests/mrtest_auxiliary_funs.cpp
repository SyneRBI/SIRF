/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
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
\ingroup MR
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

void sirf::set_unit_dcf(MRAcquisitionData& ad)
{
    std::vector<float> dcw(ad.number());
    std::fill(dcw.begin(), dcw.end(), 1.f);
    ad.set_user_floats(&dcw[0], 0);
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

void sirf::write_imagevector_to_raw(const std::string& fname_prefix, const sirf::GadgetronImagesVector& iv)
{
    for(int i=0; i<iv.items(); ++i)
    {
        std::stringstream fname_output_img;
        fname_output_img << "output_" << fname_prefix << "_image_" << i;
        write_cfimage_to_raw(fname_output_img.str(), iv.image_wrap(i));
    }
}



void sirf::set_acq_default_orientation(std::string path_in, std::string path_out)
{
    std::shared_ptr<MRAcquisitionData> sptr_ad(new AcquisitionsVector);
    AcquisitionsVector& av = (AcquisitionsVector&)*sptr_ad;
    av.read(path_in);
    int na = av.number();
    int acq_dim[10];
    av.get_acquisitions_dimensions((size_t)acq_dim);

    ISMRMRD::Acquisition acq;
    for (int i = 0; i < na; i++) {
        av.get_acquisition(i, acq);
        float* read_dir = acq.read_dir();
        float* phase_dir = acq.phase_dir();
        float* slice_dir = acq.slice_dir();

        read_dir[0] = 1.0f;
        read_dir[1] = 0.0f;
        read_dir[2] = 0.0f;
        phase_dir[0] = 0.0f;
        phase_dir[1] = 1.0f;
        phase_dir[2] = 0.0f;
        slice_dir[0] = 0.0f;
        slice_dir[1] = 0.0f;
        slice_dir[2] = 1.0f;

        av.set_acquisition(i, acq);
    }
    av.write(path_out);
}


sirf::MRAcquisitionModel
sirf::get_prepared_MRAcquisitionModel(const MRAcquisitionData& ad)
{
    sirf::GadgetronImagesVector iv;
    std::shared_ptr<GadgetronImageData> sptr_iv = std::move(iv.clone());

    std::shared_ptr<MRAcquisitionData> sptr_ad = std::move(ad.clone());

    sirf::CoilSensitivitiesVector csm;
    auto sptr_csm = std::make_shared<CoilSensitivitiesVector>(csm);
    sptr_csm->calculate(ad);

    // setup the acquisition model                 
    sirf::MRAcquisitionModel AM;
    AM.set_up(sptr_ad, sptr_iv);
    AM.set_csm(sptr_csm);

    return AM;
}
