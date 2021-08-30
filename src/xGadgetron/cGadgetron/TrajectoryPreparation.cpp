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

TrajPrep2D::TrajPointSet sirf::CartesianTrajectoryPrep::get_trajectory(const sirf::MRAcquisitionData& ac)
{
    if(ac.get_trajectory_type() != ISMRMRD::TrajectoryType::CARTESIAN)
        throw std::runtime_error("Please only ask to get the trajectory for acquisition data with a Cartesian trajectory.");

    if(ac.number() <= 0)
        throw std::runtime_error("Please pass a non-empty container.");

    ISMRMRD::Acquisition acq;

    TrajPrep2D::TrajPointSet traj;

    for(int ia=0; ia<ac.number(); ++ia)
    {
        ac.get_acquisition(ia, acq);

        TrajPrep2D::TrajPointType curr_point{(float)acq.idx().kspace_encode_step_1, 
                                             (float)acq.idx().kspace_encode_step_2};
        traj.push_back(curr_point);
    }

    return traj;
}


void sirf::GRPETrajectoryPrep::append_to_trajectory(TrajPointSet& tps, ISMRMRD::Acquisition& acq)
{
    if( acq.trajectory_dimensions() != 3)
        throw std::runtime_error("Please give Acquisition with a 3D RPE trajectory if you want to use it here.");

    TrajPrep3D::TrajPointType curr_point{0.f, acq.traj(1, 0), acq.traj(2, 0)}; // append only the 0th sample since the readout is cartesian for RPE
    tps.push_back(curr_point);
}

GRPETrajectoryPrep::TrajPointSet sirf::GRPETrajectoryPrep::calculate_trajectory(Acquisition& acq) const
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

    TrajPrep3D::TrajPointSet traj;

    for(size_t i_sample=0; i_sample<acq.number_of_samples();++i_sample)
    {
        TrajPrep3D::TrajPointType pt{0,
                                     pe_radius * cos(pe_angle),
                                     pe_radius * sin(pe_angle)};
        
        traj.push_back(pt);
    }

    return traj;
}

void sirf::NonCartesian2DTrajPrep::append_to_trajectory(TrajPointSet& tps, ISMRMRD::Acquisition& acq)
{
    if(acq.trajectory_dimensions() != 2)
        throw std::runtime_error("Please give Acquisition with a 3D RPE trajectory if you want to use it here.");

    for(int ns=0; ns<acq.number_of_samples(); ++ns)
    {
        TrajPrep2D::TrajPointType curr_point{acq.traj(0,ns), acq.traj(1,ns)};
        tps.push_back(curr_point);
    }
}

Radial2DTrajprep::TrajPointSet sirf::Radial2DTrajprep::calculate_trajectory(Acquisition& acq) const
{
    ISMRMRD::Limit rad_lims(0,0,0), ang_lims(0,0,0);
    
    if(this->kspace_encoding_.encodingLimits.kspace_encoding_step_0.is_present())
        rad_lims = this->kspace_encoding_.encodingLimits.kspace_encoding_step_0.get();

    if(this->kspace_encoding_.encodingLimits.kspace_encoding_step_1.is_present())
        ang_lims = this->kspace_encoding_.encodingLimits.kspace_encoding_step_1.get();
    
    const ISMRMRD::EncodingCounters idx = acq.idx();

    unsigned short num_angles = ang_lims.maximum;
    float const pe_angle = SIRF_PI/(float)num_angles * idx.kspace_encode_step_1;

    float const traj_norm = 2*std::max<float>( rad_lims.center-rad_lims.minimum, rad_lims.maximum-rad_lims.center);

    TrajPointSet traj;

    for(size_t i_sample=0; i_sample<acq.number_of_samples();++i_sample)
    {
        float pe_radius = float(i_sample - rad_lims.center);
        pe_radius /= traj_norm;

        TrajPointType pt{pe_radius * cos(pe_angle),
                         pe_radius * sin(pe_angle)};
        
        traj.push_back(pt);
    }

    return traj;    
}