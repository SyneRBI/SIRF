#ifndef ENCODING_H
#define ENCODING_H

#include <vector>

#include <ismrmrd/xml.h>
#include "sirf/Gadgetron/gadgetron_data_containers.h"
//#include "sirf/Gadgetron/gadgetron_image_wrap.h"
#include "sirf/iUtilities/LocalisedException.h"

#include <ismrmrd/ismrmrd.h>

#include <gadgetron/hoNDArray.h>
#include <gadgetron/vector_td.h>
#include <gadgetron/vector_td_utilities.h>

#include <gadgetron/hoNDFFT.h>
#include <gadgetron/hoNFFT.h>

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
\ingroup CartesianTrajectoryPrep
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
\ingroup GRPETrajectoryPrep
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
    std::vector< uint16_t > const rad_shift_ = {0, 2, 1, 3}; //this is bit-reversed {0 1 2 3}
    uint16_t circ_mod(uint16_t const a, uint16_t const b){ return (((a%b) + b ) % b);}
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
\ingroup CartesianFourierEncoding
\brief FFT between cartesian spaces
*/

class CartesianFourierEncoding : public FourierEncoding
{
public:
    CartesianFourierEncoding() : FourierEncoding() {}

    virtual void forward(MRAcquisitionData& ac, const CFImage* ptr_img);
    virtual void backward(CFImage* ptr_img, const MRAcquisitionData& ac);

};

typedef Gadgetron::hoNDArray<Gadgetron::floatd2> SirfTrajectoryType2D;

/*!
\ingroup CartesianFourierEncoding
\brief Radial phase encoding FFT
*/

class RPEFourierEncoding : public FourierEncoding
{
public:
    RPEFourierEncoding(): FourierEncoding() {}

    virtual void forward(MRAcquisitionData& ac, const CFImage* ptr_img);
    virtual void backward(CFImage* ptr_img, const MRAcquisitionData& ac);

    SirfTrajectoryType2D get_trajectory(const MRAcquisitionData& ac) const;

};


using namespace Gadgetron;

typedef Gadgetron::hoNDArray<std::complex<float> > CFGThoNDArr;
class Gridder_2D
{

public:

    Gridder_2D(std::vector<size_t> img_dims_output, const SirfTrajectoryType2D &traj) : nufft_operator_(from_std_vector<size_t, 2>(img_dims_output), (float)this->oversampling_factor_, (float)this->kernel_size_)

    {
        setup_nufft(img_dims_output, traj);
    }

    void setup_nufft(std::vector<size_t> img_dims_output, const SirfTrajectoryType2D &traj);


    void fft(CFGThoNDArr& kdata, const CFGThoNDArr& img);
    void ifft(CFGThoNDArr& img, const CFGThoNDArr& kdata);

protected:
    static const size_t oversampling_factor_ = 2;
    static size_t const kernel_size_ = 2;

    std::vector<size_t> trajdims_;
    std::vector<size_t> output_dims_;

    Gadgetron::hoNFFT_plan<float, 2> nufft_operator_;
};

} // namespace sirf
#endif // ENCODING_H
