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
\ingroup Gadgetron Extensions
\brief File for cartesian fourier encoding and trajectory setting.

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
\ingroup Gadgetron Extensions
\brief Specification file for preparing MRAcquisitionData for Fourier encoding.

\author Johannes Mayer
\author CCP PETMR
*/


namespace sirf{

/*!
\ingroup Gadgetron Extensions
\brief Abstract class defining the interface to set trajectories

* The strategy is to perform the trajectory computation as a pre-processsing step to the
* reconstruction. The ISMRMRD format has a 3D trajectory data field in their ISRMRMRD::Acquisition classe.
* The interface provides set_trajectory() which populates this data field depending on the implementation.
*/
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

/*!
\ingroup Gadgetron Extensions
\brief Class to set the gold-angle radial phase encoding (RPE) trajecotry
*
* Computation is based on doi:10.1002/mrm.22102
* The data reconstructed with this trajectory are parallel readouts arranged on
* a non-cartesian grid in the phase-slice-encoding plane.
* As a cartesian FFT is performed along the readout-direction the 0-dimension
* of the trajectory is set to 0 for all data points.
*/

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
\ingroup Gadgetron Extensions
\brief Abstract class defining the interface to perform Fourier transforms
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
\ingroup Gadgetron Extensions
\brief Class to perform a cartesian FFT

* The transform performs sorting of the data into a 4D matrix
* with the dimensions x, y, z and coil with a subsequent FFT along
* the first three dimensions.
* The methods assume that the acquisitions belong to the same dynamics,
* i.e. to the same echo, repetition etc.
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
