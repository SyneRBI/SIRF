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


namespace sirf{


typedef Gadgetron::hoNDArray<Gadgetron::floatd2> GadgetronTrajectoryType2D;

/*!
\ingroup Gadgetron Extensions
\brief Implementation to perform a non-cartesian FFT for RPE MR data
*
* backward(...): first data are sorted into a 4D matrix consisting of kx, ky, kz, coil
* and then the Fourier-transformed along the kx dimension. In the second step a NUFFT
* along the two remaining dimensions is performed.
*
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
\ingroup Gadgetron Extensions
\brief Class to perform a NUFFT for 2D data
*
*
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
