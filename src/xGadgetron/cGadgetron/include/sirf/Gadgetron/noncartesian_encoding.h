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

This file contains code snippets from
Gadgetron/apps/clients/gadgetron_ismrmrd_client/gadgetron_ismrmrd_client.cpp
by Michael S. Hansen

GADGETRON SOFTWARE LICENSE V1.0, NOVEMBER 2011

PERMISSION IS HEREBY GRANTED, FREE OF CHARGE, TO ANY PERSON OBTAINING
A COPY OF THIS SOFTWARE AND ASSOCIATED DOCUMENTATION FILES (THE
"SOFTWARE"), TO DEAL IN THE SOFTWARE WITHOUT RESTRICTION, INCLUDING
WITHOUT LIMITATION THE RIGHTS TO USE, COPY, MODIFY, MERGE, PUBLISH,
DISTRIBUTE, SUBLICENSE, AND/OR SELL COPIES OF THE SOFTWARE, AND TO
PERMIT PERSONS TO WHOM THE SOFTWARE IS FURNISHED TO DO SO, SUBJECT TO
THE FOLLOWING CONDITIONS:

THE ABOVE COPYRIGHT NOTICE, THIS PERMISSION NOTICE, AND THE LIMITATION
OF LIABILITY BELOW SHALL BE INCLUDED IN ALL COPIES OR REDISTRIBUTIONS
OF SUBSTANTIAL PORTIONS OF THE SOFTWARE.

SOFTWARE IS BEING DEVELOPED IN PART AT THE NATIONAL HEART, LUNG, AND BLOOD
INSTITUTE, NATIONAL INSTITUTES OF HEALTH BY AN EMPLOYEE OF THE FEDERAL
GOVERNMENT IN THE COURSE OF HIS OFFICIAL DUTIES. PURSUANT TO TITLE 17,
SECTION 105 OF THE UNITED STATES CODE, THIS SOFTWARE IS NOT SUBJECT TO
COPYRIGHT PROTECTION AND IS IN THE PUBLIC DOMAIN. EXCEPT AS CONTAINED IN
THIS NOTICE, THE NAME OF THE AUTHORS, THE NATIONAL HEART, LUNG, AND BLOOD
INSTITUTE (NHLBI), OR THE NATIONAL INSTITUTES OF HEALTH (NIH) MAY NOT
BE USED TO ENDORSE OR PROMOTE PRODUCTS DERIVED FROM THIS SOFTWARE WITHOUT
SPECIFIC PRIOR WRITTEN PERMISSION FROM THE NHLBI OR THE NIH.THE SOFTWARE IS
PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*!
\file
\ingroup SIRF non-cartesian encoding
\brief Specification file for non-cartesian Fourier encoding.

\author Johannes Mayer
*/


#ifndef NONCARTESIAN_ENCODING_H
#define NONCARTESIAN_ENCODING_H

#include <sirf/Gadgetron/FourierEncoding.h>

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
* and then Fourier-transformed along the kx dimension. In the second step a NUFFT
* along the two remaining dimensions is performed.
*
*/


class RPEFourierEncoding : public FourierEncoding
{
public:
    RPEFourierEncoding(): FourierEncoding() {}

    virtual void forward(MRAcquisitionData& ac, const CFImage& img) const;
    virtual void backward(CFImage& img, const MRAcquisitionData& ac) const;
protected:
    GadgetronTrajectoryType2D get_trajectory(const MRAcquisitionData& ac) const;

};

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

    Gridder_2D(const std::vector<size_t> img_output_dims, const GadgetronTrajectoryType2D &traj) : 
    nufft_operator_(Gadgetron::from_std_vector<size_t, 2>(img_output_dims), (float)this->oversampling_factor_, (float)this->kernel_size_)
    {
        setup_nufft(img_output_dims, traj);
    }

    void setup_nufft(const std::vector<size_t> img_output_dims, const GadgetronTrajectoryType2D &traj);


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
