#include "sirf/Gadgetron/fourierencoding.h"

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

void sirf::CartesianFourierEncoding::forward(MRAcquisitionData& ac, CFImage& img) const
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

    for (unsigned int c = 0; c < nc; c++) {
        for (unsigned int z = 0; z < nz; z++) {
            for (unsigned int y = 0; y < ny; y++) {
                for (unsigned int x = 0; x < nx; x++) {
                    ci(x, y, z, c) = (complex_float_t)img.operator()(x, y, z, c);
                }
            }
        }
    }

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

    std::vector<size_t> dims, dims_dcf;

    dims.push_back(readout);
    dims.push_back(ny); dims_dcf.push_back(ny);
    dims.push_back(nz); dims_dcf.push_back(nz);
    dims.push_back(nc);

    ISMRMRD::Limit ky_lim, kz_lim(0,0,0);

    ky_lim = e.encodingLimits.kspace_encoding_step_1.get();
    if(e.encodingLimits.kspace_encoding_step_2.is_present())
        kz_lim = e.encodingLimits.kspace_encoding_step_2.get();

    ISMRMRD::NDArray<complex_float_t> ci(dims);
    ISMRMRD::NDArray<float> dcf(dims_dcf);

    memset(ci.getDataPtr(), 0, ci.getDataSize());
    memset(dcf.getDataPtr(), 0, dcf.getDataSize());

    for (int a=0; a < ac.number(); a++) {
        ac.get_acquisition(a, acq);
        int y = ny/2 - ky_lim.center + acq.idx().kspace_encode_step_1 ;
        int z = nz/2 - kz_lim.center + acq.idx().kspace_encode_step_2;
        dcf(y, z) += (float)1;
        for (unsigned int c = 0; c < nc; c++) {
            for (unsigned int s = 0; s < readout; s++) {
                ci(s, y, z, c) += acq.data(s, c);
            }
        }
    }

    // correct for multi-acquired PE points
    for(unsigned int c=0; c<nc; ++c)
    for(unsigned int z=0; z<nz; ++z)
    for(unsigned int y=0; y<ny; ++y)
    for(unsigned int x=0; x<readout; ++x)
        ci(x,y,z,c) /= (complex_float_t)std::max(1.f, dcf(y,z));


    // now if image and kspace have different dimension then you need to interpolate or pad with zeros here
    ifft3c(ci);

    unsigned int ny_img = e.reconSpace.matrixSize.y;
    unsigned int nz_img = e.reconSpace.matrixSize.z;

    if( ny!=ny_img || nz!=nz_img)
        throw LocalisedException("Phase and slice encoding are not consistent between reconstructed image and k-space.", __FILE__, __LINE__);

    img.resize(nx_img, ny_img, nz_img, nc);
    memcpy((&img)->begin(), ci.begin(), ci.getDataSize());

    // set the header correctly of the image
    this->match_img_header_to_acquisition(img, acq);

}


