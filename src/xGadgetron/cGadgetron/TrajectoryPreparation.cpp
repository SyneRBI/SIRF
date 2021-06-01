/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 - 2021 Physikalisch-Technische Bundesanstalt (PTB)

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
\ingroup Gadgetron Extensions
\brief Utilities for setting trajectories in ISMRMRD data.
\author Johannes Mayer

*/

#include "sirf/Gadgetron/TrajectoryPreparation.h"

#include "sirf/iUtilities/LocalisedException.h"

#include <sstream>
#include <math.h>

const static double SIRF_PI = 3.14159265358979323846;
const static double SIRF_GOLDEN_ANGLE = SIRF_PI*0.618034;

using namespace sirf;
using namespace ISMRMRD;

void sirf::TrajectoryPreparation::update_acquisitions_info(MRAcquisitionData& mr_acq)
{

    IsmrmrdHeader hdr = mr_acq.acquisitions_info().get_IsmrmrdHeader();

    if(hdr.encoding.size() != 1)
        throw LocalisedException("Currrently only files with one encoding are supported", __FILE__, __LINE__);

    hdr.encoding[0].trajectory = this->traj_type_;

    this->kspace_encoding_ = hdr.encoding[0];

    std::stringstream hdr_stream;
    serialize(hdr, hdr_stream);

    AcquisitionsInfo ai(hdr_stream.str());
    mr_acq.set_acquisitions_info(ai);
}

SIRFTrajectoryType2D sirf::CartesianTrajectoryPrep::get_trajectory(const sirf::MRAcquisitionData& ac)
{
    if(ac.get_trajectory_type() != ISMRMRD::TrajectoryType::CARTESIAN)
        throw std::runtime_error("Please only ask to get the trajectory for acquisition data with a Cartesian trajectory.");

    if(ac.number() <= 0)
        throw std::runtime_error("Please pass a non-empty container.");

    ISMRMRD::Acquisition acq;

    SIRFTrajectoryType2D traj;

    for(int ia=0; ia<ac.number(); ++ia)
    {
        ac.get_acquisition(ia, acq);

        std::pair<float, float> curr_point;
        curr_point.first = acq.idx().kspace_encode_step_1;
        curr_point.second = acq.idx().kspace_encode_step_2;

        traj.push_back(curr_point);
    }

    return traj;
}


void sirf::GRPETrajectoryPrep::set_trajectory(MRAcquisitionData& mr_acq)
{
    update_acquisitions_info(mr_acq);

    for(size_t ia=0; ia<mr_acq.number(); ++ia)
    {
        Acquisition acq;
        mr_acq.get_acquisition(ia, acq);
        this->set_acquisition_trajectory(acq);
        mr_acq.set_acquisition(ia, acq);
    }
}

SIRFTrajectoryType2D sirf::GRPETrajectoryPrep::get_trajectory(const sirf::MRAcquisitionData& ac)
{
    if(ac.get_trajectory_type() != ISMRMRD::TrajectoryType::OTHER)
        throw std::runtime_error("Please only ask to get the trajectory for acquisition data with an RPE trajectory pre-computed in the acquisitions.");

    if(ac.number() <= 0)
        throw std::runtime_error("Please pass a non-empty container.");

    ISMRMRD::Acquisition acq;
    ac.get_acquisition(0, acq);

    if( acq.trajectory_dimensions() != 3)
        throw std::runtime_error("Please give Acquisition with a 3D RPE trajectory if you want to use it here.");

    SIRFTrajectoryType2D traj;

    for(int ia=0; ia<ac.number(); ++ia)
    {
        ac.get_acquisition(ia, acq);

        std::pair<float, float> curr_point;
        curr_point.first = acq.traj(1, 0);
        curr_point.second = acq.traj(2, 0);

        traj.push_back(curr_point);
    }

    return traj;
}

void sirf::GRPETrajectoryPrep::set_acquisition_trajectory(Acquisition& acq) const
{
    acq.resize(acq.number_of_samples(),acq.active_channels(), this->traj_dim_);
    std::vector<float> acq_traj = this->calculate_trajectory(acq);
    acq.setTraj(&acq_traj[0]);
}

std::vector<float> sirf::GRPETrajectoryPrep::calculate_trajectory(Acquisition& acq) const
{
    ISMRMRD::Limit rad_lims(0,0,0), ang_lims(0,0,0);
    if(this->kspace_encoding_.encodingLimits.kspace_encoding_step_1.is_present())
        rad_lims = this->kspace_encoding_.encodingLimits.kspace_encoding_step_1.get();
    if(this->kspace_encoding_.encodingLimits.kspace_encoding_step_2.is_present())
        ang_lims = this->kspace_encoding_.encodingLimits.kspace_encoding_step_2.get();

    const ISMRMRD::EncodingCounters idx = acq.idx();

//    float const pe_angle = SIRF_GOLDEN_ANGLE * (idx.kspace_encode_step_2 - 1); // for old in vivo-data
    float const pe_angle = SIRF_GOLDEN_ANGLE * idx.kspace_encode_step_2;

    size_t const num_diff_shifts = this->rad_shift_.size();
    float const rad_shift = float( this->rad_shift_.at( this->circ_mod(idx.kspace_encode_step_2 - ang_lims.center,num_diff_shifts)) ) / float(num_diff_shifts);

    float pe_radius = idx.kspace_encode_step_1 - rad_lims.center;
    pe_radius = (pe_radius==0) ? pe_radius : pe_radius+rad_shift;

    float const traj_norm = 2*std::max<float>(( rad_lims.center - rad_lims.minimum + 0), (rad_lims.maximum - rad_lims.center + (num_diff_shifts-1)/num_diff_shifts));
    pe_radius /= traj_norm;

    std::vector<float> traj;

    for(size_t i_sample=0; i_sample<acq.number_of_samples();++i_sample)
    {
        traj.push_back(0); //dummy for RPE as the readout is cartesian
        traj.push_back(pe_radius * cos( pe_angle ));
        traj.push_back(pe_radius * sin( pe_angle ));
    }

    return traj;
}
