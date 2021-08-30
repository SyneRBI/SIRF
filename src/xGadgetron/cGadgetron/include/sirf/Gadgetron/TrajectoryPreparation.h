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
template <uint16_t D>
class TrajectoryPreparation{

public:
    typedef typename std::array<float, D> TrajPointType;
    typedef typename std::vector<TrajPointType> TrajPointSet;

    TrajectoryPreparation(){}
    
    virtual void set_trajectory(MRAcquisitionData& mr_acq)
    {
        update_acquisitions_info(mr_acq);

        for(size_t ia=0; ia<mr_acq.number(); ++ia)
        {
            ISMRMRD::Acquisition acq;
            mr_acq.get_acquisition(ia, acq);
            this->set_acquisition_trajectory(acq);
            mr_acq.set_acquisition(ia, acq);
        }
    }        
    virtual TrajPointSet get_trajectory(const MRAcquisitionData& mr_acq)
    {
        if(mr_acq.number() <= 0){
            throw std::runtime_error("Please pass a non-empty container.");
        }

        ISMRMRD::Acquisition acq;
        mr_acq.get_acquisition(0, acq);

        if( acq.trajectory_dimensions() != D){
            throw std::runtime_error("Please give Acquisition with a the correct dimensionality if you want to use it here.");
        }

        TrajPointSet traj;

        for(int ia=0; ia<mr_acq.number(); ++ia)
        {
            mr_acq.get_acquisition(ia, acq);
            append_to_trajectory(traj, acq);
        }

        return traj;
    }


protected:
    ISMRMRD::Encoding kspace_encoding_;
    ISMRMRD::TrajectoryType traj_type_;

    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq) const
    {
        acq.resize(acq.number_of_samples(),acq.active_channels(), D);
        TrajPointSet acq_traj = this->calculate_trajectory(acq);
        acq.setTraj(&(acq_traj[0][0]));
    }
    virtual void update_acquisitions_info(sirf::MRAcquisitionData& mr_acq)
    {
        ISMRMRD::IsmrmrdHeader hdr = mr_acq.acquisitions_info().get_IsmrmrdHeader();

        if(hdr.encoding.size() != 1)
            throw LocalisedException("Currrently only files with one encoding are supported", __FILE__, __LINE__);

        hdr.encoding[0].trajectory = traj_type_;
        kspace_encoding_ = hdr.encoding[0];

        std::stringstream hdr_stream;
        serialize(hdr, hdr_stream);

        AcquisitionsInfo ai(hdr_stream.str());
        mr_acq.set_acquisitions_info(ai);
    }

    virtual TrajPointSet calculate_trajectory(ISMRMRD::Acquisition& acq) const =0;
    virtual void append_to_trajectory(TrajPointSet& tps, ISMRMRD::Acquisition& acq)=0;
};

typedef TrajectoryPreparation<2> TrajPrep2D;
typedef TrajectoryPreparation<3> TrajPrep3D;


/*!
\ingroup Gadgetron Extensions
\brief Class to get cartesian encoding phase encoding locations
*
* Since no computation is required for cartesian trajectories this is not
* inherited and only a static getter is made available.
*/

class CartesianTrajectoryPrep{ 
public:
    static TrajPrep2D::TrajPointSet get_trajectory(const sirf::MRAcquisitionData& ac);
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

typedef TrajectoryPreparation<3> TrajPrep3D;

class GRPETrajectoryPrep : public TrajPrep3D {

public:
    GRPETrajectoryPrep() : TrajPrep3D() {
        traj_type_ = ISMRMRD::TrajectoryType::OTHER;
    }

protected:
    TrajPointSet calculate_trajectory(ISMRMRD::Acquisition& acq) const;
    void append_to_trajectory(TrajPointSet& tps, ISMRMRD::Acquisition& acq);

private:    
    uint16_t circ_mod(uint16_t const a, uint16_t const b) const { return (((a%b) + b ) % b);}
    const std::vector< uint16_t > rad_shift_ = {0, 2, 1, 3}; //this is bit-reversed {0 1 2 3}


};
}
#endif// TRAJECTORYPREPARATION_H