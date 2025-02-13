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

#include <cstring>

#include "sirf/Gadgetron/NonCartesianEncoding.h"
#include "sirf/Gadgetron/TrajectoryPreparation.h"

using namespace sirf;
using namespace ISMRMRD;

Gridder2D::TrajectoryArrayType RPEFourierEncoding::get_trajectory(const MRAcquisitionData& ac) const
{
    sirf::GRPETrajectoryPrep tp;
    TrajectoryPreparation3D::TrajPointSet sirftraj = tp.get_trajectory(ac);

    Gridder2D::TrajectoryArrayType traj(sirftraj.size());
    traj.fill(Gadgetron::floatd2(0.f, 0.f));

    for(int ik=0; ik<traj.get_number_of_elements(); ++ik)
    {
        traj.at(ik)[0] = sirftraj[ik][1];
        traj.at(ik)[1] = sirftraj[ik][2];
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

    Gridder2D::TrajectoryArrayType traj = this->get_trajectory(ac);

    Gridder2D nufft(img_slice_dims, traj);

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

    Gridder2D::TrajectoryArrayType traj = this->get_trajectory(ac);
    size_t const num_kdata_pts = traj.get_number_of_elements();

    std::vector < size_t > img_slice_dims{img_dims[1], img_dims[2]};
    Gridder2D nufft(img_slice_dims, traj);

    std::vector< size_t> output_dims{img_dims[0], num_kdata_pts, img_dims[3]};
    CFGThoNDArr kdata(output_dims);

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

Gridder2D::TrajectoryArrayType NonCartesian2DEncoding::get_trajectory(const MRAcquisitionData& ac) const
{
    sirf::Radial2DTrajprep tp;
    TrajectoryPreparation2D::TrajPointSet sirftraj = tp.get_trajectory(ac);

    Gridder2D::TrajectoryArrayType traj(sirftraj.size());
    traj.fill(Gadgetron::floatd2(0.f, 0.f));

    for(int ik=0; ik<traj.get_number_of_elements(); ++ik)
    {
        traj.at(ik)[0] = sirftraj[ik][0];
        traj.at(ik)[1] = sirftraj[ik][1];
    }

    return traj;
}

std::vector<int> NonCartesian2DEncoding::get_slice_encoding_subset_indices(const MRAcquisitionData& full_dataset, unsigned int kspace_enc_step_2) const
{
    ISMRMRD::Acquisition acq;
    std::vector<int> subset_indices;
    for(int i=0; i<full_dataset.number(); ++i)
    {   
        full_dataset.get_acquisition(i,acq);
        if ((unsigned int)acq.idx().kspace_encode_step_2 == kspace_enc_step_2)
            subset_indices.push_back(i);
    }
    
    return subset_indices;
}

void NonCartesian2DEncoding::forward(MRAcquisitionData& ac, const CFImage& img) const 
{
    ASSERT( ac.number() >0, "Give a non-empty rawdata container if you want to use forward.");
    
    ISMRMRD::IsmrmrdHeader hdr = ac.acquisitions_info().get_IsmrmrdHeader();
    EncodingSpace enc_space = hdr.encoding[0].encodedSpace;

    ASSERT(img.getMatrixSizeZ() == enc_space.matrixSize.z, 
           "The number of slices in encoded space and image differ. Please give a slice-consistent rawdata file.");


    ISMRMRD::TrajectoryType traj_in_rawdata = ac.get_trajectory_type();
    ASSERT(traj_in_rawdata == ISMRMRD::TrajectoryType::RADIAL || 
           traj_in_rawdata == ISMRMRD::TrajectoryType::GOLDENANGLE ||
           traj_in_rawdata == ISMRMRD::TrajectoryType::SPIRAL , 
           "Give a MRAcquisitionData reference with the trajectory type RADIAL, GOLDENANGLE or SPIRAL.");

    uint16_t  Nx = img.getMatrixSizeX();
    uint16_t  Ny = img.getMatrixSizeY();
    uint16_t  NSlice = img.getMatrixSizeZ();
    uint16_t  NChannel = img.getNumberOfChannels();

    std::vector<size_t> img_dims{NSlice, Nx, Ny, NChannel};
    CFGThoNDArr img_data(img_dims);
    
    for(int nc=0; nc<NChannel; nc++)
    for(int ns=0; ns<NSlice; ns++)
    for(int ny=0; ny<Ny; ny++)
    for(int nx=0; nx<Nx; nx++)
    {   int const access_idx =  nx + Nx * (ny + Ny*(ns + nc*NSlice));
        img_data(ns,nx,ny,nc) = *(img.getDataPtr() + access_idx);
    }
    float const fft_normalisation_factor = sqrt((float)NSlice);
    Gadgetron::hoNDFFT< float >::instance()->fft1c(img_data);

    std::vector < size_t > img_slice_dims{Nx, Ny};

//    #pragma omp parallel
    for(size_t islice=0; islice < NSlice; ++islice)
    {
        std::vector<int> index_acqs_for_this_slice = ac.get_slice_encoding_index(islice);            
        std::unique_ptr<MRAcquisitionData> uptr_slice_subset = ac.clone();
        uptr_slice_subset->empty();
        ac.get_subset(*uptr_slice_subset, index_acqs_for_this_slice);

        Gridder2D::TrajectoryArrayType traj = this->get_trajectory(*uptr_slice_subset);
        Gridder2D nufft(img_slice_dims, traj);
        const size_t num_kdata_pts = traj.get_number_of_elements();

        const std::vector< size_t> output_dims{num_kdata_pts,NChannel};
        CFGThoNDArr kdata(output_dims);

        for(size_t ichannel=0; ichannel < NChannel; ++ichannel)
        {
            CFGThoNDArr img_slice(img_slice_dims);
            
            for(int ny=0; ny < Ny; ++ny)
            for(int nx=0; nx < Nx; ++nx)
                img_slice(nx,ny)= img_data(islice,nx,ny,ichannel);

            CFGThoNDArr k_slice_data_sausage;
            nufft.fft(k_slice_data_sausage, img_slice);
            
            for( int ik=0; ik<num_kdata_pts; ++ik)
                kdata(ik, ichannel) = k_slice_data_sausage.at(ik);
        }

        ISMRMRD::Acquisition acq;

        for(int ia=0; ia<uptr_slice_subset->number(); ++ia)
        {
            uptr_slice_subset->get_acquisition(ia, acq);

            for(int is=0; is<acq.number_of_samples(); ++is)
            for(int ic=0; ic<acq.active_channels(); ++ic)
            {
                const size_t access_idx = acq.number_of_samples()*ia + is;
                acq.data(is,ic) = fft_normalisation_factor * kdata(access_idx, ic);
            }
            uptr_slice_subset->set_acquisition(ia, acq);
        }
        ac.set_subset(*uptr_slice_subset, index_acqs_for_this_slice);
    }
}

void NonCartesian2DEncoding::backward(CFImage& img, const MRAcquisitionData& ac) const
{
    ISMRMRD::TrajectoryType traj_in_rawdata = ac.get_trajectory_type();
    
    ASSERT(traj_in_rawdata == ISMRMRD::TrajectoryType::RADIAL || 
           traj_in_rawdata == ISMRMRD::TrajectoryType::GOLDENANGLE ||
           traj_in_rawdata == ISMRMRD::TrajectoryType::SPIRAL, 
           "Give a MRAcquisitionData reference with the trajectory type RADIAL or GOLDENANGLE.");

    ISMRMRD::IsmrmrdHeader hdr = ac.acquisitions_info().get_IsmrmrdHeader();
    EncodingSpace rec_space = hdr.encoding[0].reconSpace;
    EncodingSpace enc_space = hdr.encoding[0].encodedSpace;

    ASSERT(rec_space.matrixSize.z == enc_space.matrixSize.z, 
           "The number of slices in encoded and reconstructed space differ. Please give rawdata from a slice-consistent file.");

    std::vector<size_t> kspace_dims;
    ac.get_kspace_dimensions(kspace_dims);

    size_t const Nx     = rec_space.matrixSize.x ;
    size_t const Ny     = rec_space.matrixSize.y ;
    size_t const NSlice = enc_space.matrixSize.z;
    size_t const NChannel = kspace_dims[3];
    
    float const fft_normalisation_factor = sqrt(float(NSlice));

    auto uptr_slice_subset = ac.clone();

    std::vector<size_t> img_dimensions{NSlice,Nx,Ny,NChannel};
    CFGThoNDArr img_data(img_dimensions);

    for(size_t islice=0; islice<NSlice; ++islice)
    {
        uptr_slice_subset->empty();

        std::vector<int> slice_subset_indices = get_slice_encoding_subset_indices(ac, islice);
        ac.get_subset(*uptr_slice_subset, slice_subset_indices);

        Gridder2D::TrajectoryArrayType traj = this->get_trajectory(*uptr_slice_subset);
        const size_t num_kdata_pts = traj.get_number_of_elements();
        std::vector<size_t> const kdata_dims{num_kdata_pts, NChannel};
        CFGThoNDArr kspace_data(kdata_dims);
        
        for(int ia=0; ia<uptr_slice_subset->number(); ++ia)
        {
            ISMRMRD::Acquisition acq;
            uptr_slice_subset->get_acquisition(ia, acq);

            for(int nc=0; nc<acq.active_channels(); ++nc)
            for(int ns=0; ns<acq.number_of_samples(); ++ns)
            {   
                int const access_idx = kspace_dims[0] * ia + ns;
                kspace_data(access_idx, nc) = acq.data(ns,nc);
            }
        }
       
        std::vector < size_t > img_slice_dims{Nx, Ny};
        Gridder2D nufft(img_slice_dims, traj);

        for(size_t ichannel=0; ichannel<NChannel; ++ichannel)
        {
            CFGThoNDArr k_slice_data_sausage(num_kdata_pts);
            k_slice_data_sausage.fill(complex_float_t(0,0));

            for(int ik=0;ik<num_kdata_pts;++ik)
            {
                k_slice_data_sausage.at(ik) = fft_normalisation_factor*kspace_data(ik,ichannel);
            }

            CFGThoNDArr imgdata_slice;
            nufft.ifft(imgdata_slice, k_slice_data_sausage);

            for(size_t iy=0; iy<Ny; ++iy)
            for(size_t ix=0; ix<Nx; ++ix)
                img_data(islice, ix, iy, ichannel) = imgdata_slice(ix, iy);
        }
    }

    Gadgetron::hoNDFFT< float >::instance()->ifft1c(img_data);

    img.resize(Nx, Ny, NSlice, NChannel);

    for(int nc=0; nc<NChannel; ++nc)
    for(int ny=0; ny<Ny; ++ny)
    for(int nx=0; nx<Nx; ++nx)
    for(int ns=0; ns<NSlice; ++ns)
        img(nx,ny,ns,nc) = img_data(ns,nx,ny,nc);
    
    img.setFieldOfView( rec_space.fieldOfView_mm.x, rec_space.fieldOfView_mm.y ,rec_space.fieldOfView_mm.z );

    ISMRMRD::Acquisition acq;
    ac.get_acquisition(0,acq);
    this->match_img_header_to_acquisition(img, acq);    
}
