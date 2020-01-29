/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2019 University College London

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

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
\ingroup Registration
\brief NiftyReg's aladin class for rigid and affine registrations.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/SPM12Registration.h"
#include "sirf/Reg/NiftiImageData.h"
#include <sys/stat.h>
#include <MatlabEngine.hpp>
#include <boost/filesystem.hpp>

using namespace sirf;
using namespace matlab::engine;
using namespace matlab::data;

inline bool check_file_exists(const std::string& filename, const bool existance_allowed) {
    struct stat buffer;
    const bool file_exists (stat (filename.c_str(), &buffer) == 0);
    if (file_exists && !existance_allowed)
        throw std::runtime_error("SPM12Registration<dataType>::process(): file exists: " + filename);
    return file_exists;
}

template<class dataType>
void convert_to_NiftiImageData_if_not_already(std::shared_ptr<const NiftiImageData<dataType> > &output_sptr, const std::shared_ptr<const ImageData> &input_sptr)
{
    // Try to dynamic cast from ImageData to (const) NiftiImageData. This will only succeed if original type was NiftiImageData
    output_sptr = std::dynamic_pointer_cast<const NiftiImageData<dataType> >(input_sptr);
    // If output is a null pointer, it means that a different image type was supplied (e.g., STIRImageData).
    // In this case, construct a NiftiImageData
    if (!output_sptr)
        output_sptr = std::make_shared<const NiftiImageData<dataType> >(*input_sptr);
}

template<class dataType>
void SPM12Registration<dataType>::set_working_folder(const std::string &working_folder)
{
    // Make sure it's absolute
    _working_folder = boost::filesystem::absolute(working_folder).string();
}

template<class dataType>
void SPM12Registration<dataType>::process()
{
    // Check the paramters that are NOT set via the parameter file have been set.
    this->check_parameters();

    // Filenames
    const std::string ref_filename = _working_folder + "/ref.nii";
    const std::string flo_filename = _working_folder + "/flo.nii";
    check_file_exists(ref_filename, _working_folder_overwrite);
    check_file_exists(flo_filename, _working_folder_overwrite);

    // Convert images to matlab::data::StructArray with structure expected by SPM
    std::shared_ptr<const NiftiImageData<dataType> > ref_nifti_sptr, flo_nifti_sptr;
    convert_to_NiftiImageData_if_not_already(ref_nifti_sptr, this->_reference_image_sptr);
    convert_to_NiftiImageData_if_not_already(flo_nifti_sptr, this->_floating_image_sptr);

    // Save to file
    ref_nifti_sptr->write(ref_filename);
    flo_nifti_sptr->write(flo_filename);

    // Start MATLAB engine synchronously
    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
    std::cout << "Started MATLAB Engine" << std::endl;

    // Create MATLAB data array factory
    ArrayFactory factory;

    // Create cell array for filenames {'ref.nii','flo.nii'}
    const unsigned num_images = 2;
    CellArray spm_filenames = factory.createCellArray({num_images},
        factory.createCharArray(ref_filename),
        factory.createCharArray(flo_filename));

    // Create struct array for parameters: struct('quality',1,'rtm',1))
    StructArray spm_params = factory.createStructArray({1}, { "quality", "rtm" });
    spm_params[0]["quality"] = factory.createScalar<int>(1);
    spm_params[0]["rtm"] = factory.createScalar<int>(1);

    // Create a vector of input arguments
    std::vector<Array> args({
        spm_filenames,
        spm_params
    });

    // Call spm_realign
    const size_t num_returned = 0;
    matlabPtr->feval(u"spm_realign", num_returned, args);

    // Read the transformation matrix back in
    // Read text file
    std::string line;
    std::ifstream myfile(_working_folder + "/rp_flo.txt");
    if (!myfile.is_open())
        throw std::runtime_error("SPM12Registration::process() failed to open spm_realign results here: " + _working_folder + "/rp_flo.txt");
    try {
        getline (myfile,line);
    }
    catch (...) {
        throw std::runtime_error("SPM12Registration::process() failed to read spm_realign results here: " + _working_folder + "/rp_flo.txt");
    }
    myfile.close();

    // Convert text to numbers
    std::stringstream ss;
    ss << line;
    std::string temp;
    float found;
    std::vector<float> results;
    while (!ss.eof()) {

        /* extracting word by word from stream */
        ss >> temp;

        /* Checking the given word is integer or not */
        if (std::stringstream(temp) >> found)
            results.push_back(found);

        /* To save from space at the end of string */
        temp = "";
    }
    for (unsigned i=0; i<results.size(); ++i)
        std::cout << "restult " << i << ": " << results[i] << "\n";
}

template<class dataType>
void SPM12Registration<dataType>::check_parameters() const
{
    // Call base class
    Registration<dataType>::check_parameters();

    if (_working_folder.empty())
        throw std::runtime_error("SPM12Registration<dataType>::check_parameters(): Missing working folder.");
}

namespace sirf {
template class SPM12Registration<float>;
}
