/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2021 Rutherford Appleton Laboratory STFC

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
\ingroup Fourier encoding
\brief File for cartesian fourier encoding.

\author Johannes Mayer
\author SyneRBI
*/

#ifndef ENCODING_H
#define ENCODING_H

#include <vector>

#include <ismrmrd/xml.h>
#include "sirf/Gadgetron/gadgetron_data_containers.h"
//#include "sirf/Gadgetron/gadgetron_image_wrap.h"
#include "sirf/iUtilities/LocalisedException.h"

#include <ismrmrd/ismrmrd.h>

/*!
\file
\ingroup Fourier Encoding
\brief Specification file for preparing MRAcquisitionData for Fourier encoding.

\author Johannes Mayer
\author CCP PETMR
*/

/*!
\ingroup aTrajectoryPreparation
\brief Abstract class for trajectory preparation

*/

namespace sirf{

class aTrajectoryPreparation{

public:
    aTrajectoryPreparation(){}
    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq) =0;

protected:

    void update_acquisitions_info(sirf::MRAcquisitionData& mr_acq);

    ISMRMRD::Encoding kspace_encoding_;
    ISMRMRD::TrajectoryType traj_type_;

    uint16_t traj_dim_;
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq) const =0;
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq) const =0;
};


typedef std::vector< std::pair<float, float> > SIRFTrajectoryType2D;

class GRPETrajectoryPrep : public aTrajectoryPreparation {

public:
    GRPETrajectoryPrep(): aTrajectoryPreparation() {
        traj_type_ = ISMRMRD::TrajectoryType::OTHER;
        traj_dim_ = 3;
    }

    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq);
    static SIRFTrajectoryType2D get_trajectory(const sirf::MRAcquisitionData& mr_acq);

protected:
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq) const;
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq) const;
    uint16_t circ_mod(uint16_t const a, uint16_t const b) const { return (((a%b) + b ) % b);}

    const std::vector< uint16_t > rad_shift_ = {0, 2, 1, 3}; //this is bit-reversed {0 1 2 3}

};

/*!
\ingroup Fourier Encoding
\brief Abstract class for doing FFTs for different trajectories for self-consistent k-space data.

*/

class FourierEncoding
{
public:
    FourierEncoding(){}

    virtual void forward(MRAcquisitionData& ac, const CFImage* ptr_img)=0;
    virtual void backward(CFImage* ptr_img, const MRAcquisitionData& ac)=0;

    void match_img_header_to_acquisition(CFImage& img, const ISMRMRD::Acquisition& acq);
};

/*!
\ingroup Fourier Encoding
\brief FFT between cartesian spaces
*/

class CartesianFourierEncoding : public FourierEncoding
{
public:
    CartesianFourierEncoding() : FourierEncoding() {}

    virtual void forward(MRAcquisitionData& ac, const CFImage* ptr_img);
    virtual void backward(CFImage* ptr_img, const MRAcquisitionData& ac);

};


} // namespace sirf
#endif // ENCODING_H
