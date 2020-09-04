#ifndef ENCODING_H
#define ENCODING_H

#include <vector>

#include <ismrmrd/xml.h>
#include "sirf/Gadgetron/gadgetron_data_containers.h"



/*!
\file
\ingroup Fourier Encoding
\brief Specification file for preparing MRAcquisitionData for Fourier encoding.

\author Johannes Mayer
\author CCP PETMR
*/

/*!
\ingroup Fourier Encoding
\brief Abstract class for trajectory preparation

*/
namespace sirf{

class aTrajectoryPreparation{

public:
    aTrajectoryPreparation(){}
    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq)=0;

protected:

    void update_acquisitions_info(sirf::MRAcquisitionData& mr_acq);

    ISMRMRD::Encoding kspace_encoding_;
    ISMRMRD::TrajectoryType traj_type_;

    uint16_t traj_dim_;
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq)=0;
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq)=0;
};

/*!
\ingroup Fourier Encoding
\brief Cartesian trajectory preparation class

*/

class CartesianTrajectoryPrep : public aTrajectoryPreparation{

public:
    CartesianTrajectoryPrep(): aTrajectoryPreparation() {
        traj_type_ = ISMRMRD::TrajectoryType::CARTESIAN;
        traj_dim_ = 0;
    }

    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq);

protected:
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq){}
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq){return std::vector<float>{};}
};

/*!
\ingroup Fourier Encoding
\brief Golden Radial Phase Encoding interleaved trajectory preparation class.

*/
class GRPETrajectoryPrep : public aTrajectoryPreparation {

public:
    GRPETrajectoryPrep(): aTrajectoryPreparation() {
        traj_type_ = ISMRMRD::TrajectoryType::OTHER;
        traj_dim_ = 3;
    }

    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq);

protected:
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq);
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq);
    std::vector< uint16_t > const rad_shift_ = {0, 2, 1, 3};
    uint16_t circ_mod(uint16_t const a, uint16_t const b){ return (((a%b) + b ) % b);}
};




class FourierEncoding
{
public:
    FourierEncoding(){}

    virtual void forward(CFImage* ptr_img, MRAcquisitionData& ac)=0;
    virtual void backward(CFImage* ptr_img, MRAcquisitionData& ac)=0;

    virtual void match_img_header_to_acquisition(CFImage& img, const ISMRMRD::Acquisition& acq);
};

class CartesianFourierEncoding : public FourierEncoding
{
public:
    CartesianFourierEncoding() : FourierEncoding() {}

    virtual void forward(CFImage* ptr_img, MRAcquisitionData& ac);
    virtual void backward(CFImage* ptr_img, MRAcquisitionData& ac);

};



} // namespace sirf
#endif // ENCODING_H
