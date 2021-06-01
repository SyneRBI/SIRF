/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020-2021 Physikalisch-Technische Bundesanstalt (PTB)

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
\brief Implementation file for non-cartesian Fourier encoding.

\author Johannes Mayer
*/

#include "sirf/Gadgetron/NonCartesianEncoding.h"
#include "sirf/Gadgetron/TrajectoryPreparation.h"

using namespace sirf;
using namespace ISMRMRD;

GadgetronTrajectoryType2D RPEFourierEncoding::get_trajectory(const MRAcquisitionData& ac) const
{
    SIRFTrajectoryType2D sirftraj = GRPETrajectoryPrep::get_trajectory(ac);

    GadgetronTrajectoryType2D traj(sirftraj.size());
    traj.fill(Gadgetron::floatd2(0.f, 0.f));

    for(int ik=0; ik<traj.get_number_of_elements(); ++ik)
    {
        traj.at(ik)[0] = sirftraj[ik].first;
        traj.at(ik)[1] = sirftraj[ik].second;
    }

    return traj;
}

void RPEFourierEncoding::backward(CFImage& img, const MRAcquisitionData& ac) const
{
    ASSERT( ac.get_trajectory_type() == ISMRMRD::TrajectoryType::OTHER, "Give a MRAcquisitionData reference with the trajectory type OTHER.");

    ISMRMRD::IsmrmrdHeader hdr = ac.acquisitions_info().get_IsmrmrdHeader();
    ISMRMRD::Encoding e = hdr.encoding[0];

    std::vector<size_t> kspace_dims;
    ac.get_kspace_dimensions(kspace_dims);

    std::vector<size_t> const kdata_dims{kspace_dims[0], ac.number(), kspace_dims[3]};

    CFGThoNDArr kspace_data(kdata_dims);

    for(int ia=0; ia<ac.number(); ++ia)
    {
        ISMRMRD::Acquisition acq;
        ac.get_acquisition(ia, acq);

        for(int is=0; is<acq.number_of_samples(); ++is)
            for(int ic=0; ic<acq.active_channels(); ++ic)
                kspace_data(is, ia, ic) =  acq.data(is,ic);
    }

    Gadgetron::hoNDFFT< float >::instance()->ifft1c(kspace_data);

    EncodingSpace rec_space = e.reconSpace;
    std::vector < size_t > img_slice_dims{rec_space.matrixSize.y, rec_space.matrixSize.z};

    GadgetronTrajectoryType2D traj = this->get_trajectory(ac);

    Gridder_2D nufft(img_slice_dims, traj);

    img.resize(rec_space.matrixSize.x, rec_space.matrixSize.y, rec_space.matrixSize.z, kspace_dims[3]);

    float const fft_normalisation_factor = sqrt(float(rec_space.matrixSize.x));

    for(size_t ichannel=0; ichannel<kspace_dims[3]; ++ichannel)
    {
        for(size_t islice=0;islice<kspace_dims[0]; ++islice)
        {
            CFGThoNDArr k_slice_data_sausage(kdata_dims[1]);
            k_slice_data_sausage.fill(complex_float_t(0,0));

            for(int ik=0;ik<kdata_dims[1];++ik)
            {
                k_slice_data_sausage.at(ik) = kspace_data(islice,ik,ichannel);
            }

            CFGThoNDArr imgdata_slice;
            nufft.ifft(imgdata_slice, k_slice_data_sausage);

            for(size_t iz=0; iz<rec_space.matrixSize.z; ++iz)
            for(size_t iy=0; iy<rec_space.matrixSize.y; ++iy)
                img.operator()(islice, iy, iz, ichannel) = fft_normalisation_factor * imgdata_slice(iy, iz);
        }
    }

    img.setFieldOfView( rec_space.fieldOfView_mm.x, rec_space.fieldOfView_mm.y ,rec_space.fieldOfView_mm.z );

    ISMRMRD::Acquisition acq;
    ac.get_acquisition(0,acq);

    this->match_img_header_to_acquisition(img, acq);
}

void RPEFourierEncoding::forward(MRAcquisitionData& ac, const CFImage& img) const
{

    ASSERT( ac.number() >0, "Give a non-empty rawdata container if you want to use the rpe forward.");
    ASSERT( ac.get_trajectory_type() == ISMRMRD::TrajectoryType::OTHER, "Give a MRAcquisitionData reference with the trajectory type OTHER.");

    std::vector<size_t> img_dims;
    img_dims.push_back(img.getMatrixSizeX());
    img_dims.push_back(img.getMatrixSizeY());
    img_dims.push_back(img.getMatrixSizeZ());
    img_dims.push_back(img.getNumberOfChannels());

    CFGThoNDArr img_data(img_dims);
    std::memcpy(img_data.begin(), img.getDataPtr(), img.getDataSize());

    GadgetronTrajectoryType2D traj = this->get_trajectory(ac);
    size_t const num_kdata_pts = traj.get_number_of_elements();

    std::vector < size_t > img_slice_dims{img_dims[1], img_dims[2]};
    Gridder_2D nufft(img_slice_dims, traj);

    std::vector< size_t> output_dims{img_dims[0], num_kdata_pts, img_dims[3]};
    CFGThoNDArr kdata(output_dims);

//    #pragma omp parallel
    for(size_t ichannel=0; ichannel<img_dims[3]; ++ichannel)
    for(size_t islice=0;islice<img_dims[0]; ++islice)
    {

        CFGThoNDArr img_slice(img_slice_dims);

        for(int ny=0; ny<img_dims[1]; ++ny)
        for(int nz=0; nz<img_dims[2]; ++nz)
             img_slice(ny,nz)= img_data(islice,ny,nz,ichannel);

        CFGThoNDArr k_slice_data_sausage;
        nufft.fft(k_slice_data_sausage, img_slice);

        for( int ik=0; ik<num_kdata_pts; ++ik)
            kdata(islice, ik, ichannel) = k_slice_data_sausage.at(ik);
    }

    Gadgetron::hoNDFFT< float >::instance()->fft1c(kdata);

    ISMRMRD::Acquisition acq;
    ac.get_acquisition(0, acq);

    ASSERT( acq.number_of_samples() == img_dims[0],"NUMBER OF SAMPLES OF RAWDATA DONT MATCH IMAGES SLICES");
    ASSERT( acq.active_channels() == img_dims[3],"NUMBER OF CHANNELS OF RAWDATA DONT MATCH IMAGES CHANNELS");

    float const fft_normalisation_factor = sqrt((float)acq.number_of_samples());

    for(int ia=0; ia<num_kdata_pts; ++ia)
    {
        ac.get_acquisition(ia, acq);

        for(int is=0; is<acq.number_of_samples(); ++is)
            for(int ic=0; ic<acq.active_channels(); ++ic)
                acq.data(is,ic) = fft_normalisation_factor * kdata(is, ia, ic);

        ac.set_acquisition(ia, acq);
    }
}


void Gridder_2D::setup_nufft(const std::vector<size_t> img_output_dims, const GadgetronTrajectoryType2D &traj)
{
    if( img_output_dims.size() != 2)
        throw LocalisedException("The image dimensions of the output should be of size 2." , __FILE__, __LINE__);

    traj.get_dimensions(this->trajdims_);

    this->output_dims_ = img_output_dims;

    this->nufft_operator_.preprocess(traj);
}



void Gridder_2D::ifft(CFGThoNDArr& img, const CFGThoNDArr& kdata)
{
    auto sptr_const_dcw = std::make_shared<Gadgetron::hoNDArray<float> >( this->trajdims_);
    float const normed_dcw_value = 1.0;

    sptr_const_dcw ->fill(normed_dcw_value);

    img.create(this->output_dims_);
    img.fill(std::complex<float>(0.f, 0.f));

    this->nufft_operator_.compute(kdata, img, sptr_const_dcw.get(), Gadgetron::NFFT_comp_mode::BACKWARDS_NC2C);

}

void Gridder_2D::fft(CFGThoNDArr& kdata, const CFGThoNDArr& img)
{
    auto sptr_unit_dcw = std::make_shared<Gadgetron::hoNDArray<float> >( this->trajdims_);
    sptr_unit_dcw ->fill(1.f);

    kdata.create(this->trajdims_);

    this->nufft_operator_.compute(img, kdata, sptr_unit_dcw.get(), Gadgetron::NFFT_comp_mode::FORWARDS_C2NC);

}

