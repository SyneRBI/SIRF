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

*/

#include "sirf/Gadgetron/FourierEncoding.h"

#include <sstream>
#include <math.h>

#include "sirf/iUtilities/LocalisedException.h"

using namespace sirf;
// using namespace ISMRMRD;

void sirf::FourierEncoding::match_img_header_to_acquisition(CFImage& img, const ISMRMRD::Acquisition& acq) const
{

    auto acq_hdr = acq.getHead();
    auto idx = acq_hdr.idx;

    img.setAverage(idx.average);
    img.setSlice(idx.slice);
    img.setContrast(idx.contrast);
    img.setPhase(idx.phase);
    img.setRepetition(idx.repetition);
    img.setSet(idx.set);

    img.setReadDirection(acq_hdr.read_dir[0], acq_hdr.read_dir[1], acq_hdr.read_dir[2]);
    img.setPhaseDirection(acq_hdr.phase_dir[0], acq_hdr.phase_dir[1], acq_hdr.phase_dir[2]);
    img.setSliceDirection(acq_hdr.slice_dir[0], acq_hdr.slice_dir[1], acq_hdr.slice_dir[2]);

    img.setPosition(acq_hdr.position[0], acq_hdr.position[1], acq_hdr.position[2]);
    img.setPatientTablePosition(acq_hdr.patient_table_position[0], acq_hdr.patient_table_position[1], acq_hdr.patient_table_position[2]);

}

/*
The next two methods:

	CartesianFourierEncoding::forward
	CartesianFourierEncoding::backward

contain code snippets from ISMRMRD/utilities/generate_cartesian_shepp_logan.cpp
by Michael S. Hansen

ISMRMRD SOFTWARE LICENSE JULY 2013

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

void sirf::CartesianFourierEncoding::forward(MRAcquisitionData& ac, const CFImage& img) const
{

    std::string par;
    ISMRMRD::IsmrmrdHeader header;
    par = ac.acquisitions_info();
    ISMRMRD::deserialize(par.c_str(), header);

    if( header.encoding.size() > 1)
        LocalisedException("Currently only one encoding is supported per rawdata file.", __FUNCTION__, __LINE__);

    ISMRMRD::Encoding e = header.encoding[0];

    ISMRMRD::Acquisition acq;
    ac.get_acquisition(0, acq);


    unsigned int readout = acq.number_of_samples();
    unsigned int ny_k_space = e.encodedSpace.matrixSize.y;
    unsigned int nz_k_space = e.encodedSpace.matrixSize.z;

    ISMRMRD::Limit ky_lim, kz_lim(0,0,0);

    ky_lim = e.encodingLimits.kspace_encoding_step_1.get();
    if(e.encodingLimits.kspace_encoding_step_2.is_present())
        kz_lim = e.encodingLimits.kspace_encoding_step_2.get();



    unsigned int nx = img.getMatrixSizeX();
    unsigned int ny = img.getMatrixSizeY();
    unsigned int nz = img.getMatrixSizeZ();
    unsigned int nc = img.getNumberOfChannels();

    if(nx != readout || ny_k_space!= ny || nz_k_space != nz || nc != acq.active_channels())
        throw LocalisedException("K-space dimensions and image dimensions don't match.",   __FILE__, __LINE__);


    std::vector<size_t> dims;
    dims.push_back(nx);
    dims.push_back(ny);
    dims.push_back(nz);
    dims.push_back(nc);

    ISMRMRD::NDArray<complex_float_t> ci(dims);
    memset(ci.getDataPtr(), 0, ci.getDataSize());
    std::memcpy(ci.getDataPtr(), img.getDataPtr(), img.getDataSize());

    fft3c(ci);

    for(size_t i =0; i<ac.items(); ++i)
    {
        ac.get_acquisition(i, acq);
        acq.resize(nx, nc); // no trajectory information is set

        int ky = ny/2 - ky_lim.center + acq.idx().kspace_encode_step_1;
        int kz = nz/2 - kz_lim.center + acq.idx().kspace_encode_step_2;

        for (unsigned int c = 0; c < nc; c++) {
            for (unsigned int s = 0; s < nx; s++) {
                acq.data(s, c) = ci(s, ky, kz, c);
            }
        }
        ac.set_acquisition(i, acq);
    }
}

void sirf::CartesianFourierEncoding::backward(CFImage& img, const MRAcquisitionData& ac) const
{

    if(ac.items()<1)
        LocalisedException("No acquisitions in vector. Trying to backward transform empty vector.", __FUNCTION__, __LINE__);

    ISMRMRD::IsmrmrdHeader header = ac.acquisitions_info().get_IsmrmrdHeader();

    if( header.encoding.size() > 1)
        LocalisedException("Currently only one encoding is supported per rawdata file.", __FUNCTION__, __LINE__);

    ISMRMRD::Encoding e = header.encoding[0];

    ISMRMRD::Acquisition acq;

    ac.get_acquisition(0, acq);

    unsigned int readout = acq.number_of_samples();
    unsigned int nc = acq.active_channels();

    unsigned int ny = e.encodedSpace.matrixSize.y;
    unsigned int nz = e.encodedSpace.matrixSize.z;

    unsigned int nx_img = e.reconSpace.matrixSize.x;

    if(nx_img != readout)
        throw LocalisedException("Number of readout points and reconstructed image dimension in readout direction are assumed the same.",   __FILE__, __LINE__);

    std::vector<size_t> dims;

    dims.push_back(readout);
    dims.push_back(ny);
    dims.push_back(nz);
    dims.push_back(nc);

    ISMRMRD::Limit ky_lim, kz_lim(0,0,0);

    ky_lim = e.encodingLimits.kspace_encoding_step_1.get();
    if(e.encodingLimits.kspace_encoding_step_2.is_present())
        kz_lim = e.encodingLimits.kspace_encoding_step_2.get();

    ISMRMRD::NDArray<complex_float_t> ci(dims);
    memset(ci.getDataPtr(), 0, ci.getDataSize());
    
    for (int a=0; a < ac.number(); a++) {
        ac.get_acquisition(a, acq);
        int y = ny/2 - ky_lim.center + acq.idx().kspace_encode_step_1 ;
        int z = nz/2 - kz_lim.center + acq.idx().kspace_encode_step_2;
    
        for (unsigned int c = 0; c < nc; c++) {
            for (unsigned int s = 0; s < readout; s++) {
                ci(s, y, z, c) += acq.data(s, c);
            }
        }
    }

    // now if image and kspace have different dimension then you need to interpolate or pad with zeros here
    ifft3c(ci);

    unsigned int ny_img = e.reconSpace.matrixSize.y;
    unsigned int nz_img = e.reconSpace.matrixSize.z;

    if( ny!=ny_img || nz!=nz_img)
        throw LocalisedException("Phase and slice encoding are not consistent between reconstructed image and k-space.", __FILE__, __LINE__);

    img.resize(nx_img, ny_img, nz_img, nc);
    std::memcpy(img.getDataPtr(), ci.getDataPtr(), ci.getDataSize());

    // set the header correctly of the image
    this->match_img_header_to_acquisition(img, acq);

}


