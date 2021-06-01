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

#ifndef TRAJECTORYPREPARATION_H
#define TRAJECTORYPREPARATION_H

#include <vector>

#include <ismrmrd/xml.h>
#include <ismrmrd/ismrmrd.h>

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/iUtilities/LocalisedException.h"



namespace sirf{

/*!
\ingroup Gadgetron Extensions
\brief Abstract class defining the interface to set trajectories

* The strategy is to perform the trajectory computation as a pre-processsing step to the
* reconstruction. The ISMRMRD format has a 3D trajectory data field in their ISRMRMRD::Acquisition classe.
* The interface provides set_trajectory() which populates this data field depending on the implementation.
*/
class TrajectoryPreparation{

public:
    TrajectoryPreparation(){}
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

/*!
\ingroup Gadgetron Extensions
\brief Class to get cartesian encoding phase encoding locations
*
* Since no computation is required for cartesian trajectories this is not
* inherited and only a static getter is made available.
*/

class CartesianTrajectoryPrep{ 
public:
    static SIRFTrajectoryType2D get_trajectory(const sirf::MRAcquisitionData& ac);
};

/*!
\ingroup Gadgetron Extensions
\brief Class to set the golden-angle radial phase encoding (GRPE) trajectory
*
* Computation is based on doi:10.1002/mrm.22102
* The data reconstructed with this trajectory are parallel readouts arranged on
* a non-cartesian grid in the phase-slice-encoding plane.
* As a cartesian FFT is performed along the readout-direction the 0-dimension
* of the trajectory is set to 0 for all data points.
*/

class GRPETrajectoryPrep : public TrajectoryPreparation {

public:
    GRPETrajectoryPrep(): TrajectoryPreparation() {
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
}
#endif// TRAJECTORYPREPARATION_H