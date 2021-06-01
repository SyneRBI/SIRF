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
\brief File for cartesian fourier encoding and trajectory setting.

\author Johannes Mayer
*/

#ifndef FOURIERENCODING_H
#define FOURIERENCODING_H

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/iUtilities/LocalisedException.h"

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
\brief Abstract class defining the interface to perform Fourier transforms
*/

class FourierEncoding
{
public:
    FourierEncoding(){}

    virtual void forward(MRAcquisitionData& ac, const CFImage& img) const =0;
    virtual void backward(CFImage& img, const MRAcquisitionData& ac) const =0;
    
    void match_img_header_to_acquisition(CFImage& img, const ISMRMRD::Acquisition& acq) const;
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

    virtual void forward(MRAcquisitionData& ac, const CFImage& img) const;
    virtual void backward(CFImage& img, const MRAcquisitionData& ac) const;
    
};

} // namespace sirf
#endif // FOURIERENCODING_H