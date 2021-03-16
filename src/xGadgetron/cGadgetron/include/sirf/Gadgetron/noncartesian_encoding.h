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
\ingroup SIRF non-cartesian encoding
\brief Specification file for non-cartesian Fourier encoding.

\author Johannes Mayer
\author SyneRBI
*/


#ifndef NONCARTESIAN_ENCODING_H
#define NONCARTESIAN_ENCODING_H

#include <sirf/Gadgetron/encoding.h>

#include <gadgetron/hoNDArray.h>
#include <gadgetron/vector_td.h>
#include <gadgetron/vector_td_utilities.h>

#include <gadgetron/hoNDFFT.h>
#include <gadgetron/hoNFFT.h>

/*!
\ingroup GRPETrajectoryPrep
\brief Golden Radial Phase Encoding interleaved trajectory preparation class.
*/

namespace sirf{

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
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq);
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq);
    std::vector< uint16_t > const rad_shift_ = {0, 2, 1, 3}; //this is bit-reversed {0 1 2 3}
    uint16_t circ_mod(uint16_t const a, uint16_t const b){ return (((a%b) + b ) % b);}
};


typedef Gadgetron::hoNDArray<Gadgetron::floatd2> GadgetronTrajectoryType2D;

/*!
\ingroup RPE Fourier Encoding
\brief Radial phase encoding operator
*/

class RPEFourierEncoding : public FourierEncoding
{
public:
    RPEFourierEncoding(): FourierEncoding() {}

    virtual void forward(MRAcquisitionData& ac, const CFImage* ptr_img);
    virtual void backward(CFImage* ptr_img, const MRAcquisitionData& ac);
protected:
    GadgetronTrajectoryType2D get_trajectory(const MRAcquisitionData& ac) const;

};

using namespace Gadgetron;

typedef Gadgetron::hoNDArray<std::complex<float> > CFGThoNDArr;

/*!
\ingroup CartesianFourierEncoding
\brief Radial phase encoding operator
*/

class Gridder_2D
{
public:

    Gridder_2D(std::vector<size_t> img_dims_output, const GadgetronTrajectoryType2D &traj) : nufft_operator_(from_std_vector<size_t, 2>(img_dims_output), (float)this->oversampling_factor_, (float)this->kernel_size_)

    {
        setup_nufft(img_dims_output, traj);
    }

    void setup_nufft(std::vector<size_t> img_dims_output, const GadgetronTrajectoryType2D &traj);


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
#endif // NONCARTESIAN_ENCODING_H
