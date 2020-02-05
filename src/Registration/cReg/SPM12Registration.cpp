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
#include "sirf/Reg/AffineTransformation.h"
#include <codecvt>

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
    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB(std::vector<String>({u"-nojvm"}));
    std::cout << "Started MATLAB Engine" << std::endl;

    // We'll need to be able to convert commands from std::string to std::u16string
    std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t> convertor;

    // Call spm_realign
    const std::string spm_realign_command = "spm_output = spm_realign(char('" + ref_filename + "','" + flo_filename + "'),struct('quality',1,'rtm',0));";
    matlabPtr->eval(convertor.from_bytes(spm_realign_command));

    // Call spm_reslice
    const std::string spm_reslice_command = "spm_reslice(spm_output(2),struct('which',[2,0]));";
    matlabPtr->eval(convertor.from_bytes(spm_reslice_command));

    const std::string resliced_filename = _working_folder + "/rflo.nii";
    this->_warped_image_sptr = std::make_shared<NiftiImageData<dataType> >(resliced_filename);

    // Clean up
    if (_delete_temp_files) {
        remove( ref_filename.c_str() );
        remove( flo_filename.c_str() );
        remove( resliced_filename.c_str() );
    }

    // Get the transformation matrix
    const std::string tm_in_nr_space_command = "tm_in_nr_space = inv(spm_output(2).mat / spm_output(2).private.mat0);";
    matlabPtr->eval(convertor.from_bytes(tm_in_nr_space_command));
    TypedArray<double> tm = matlabPtr->getVariable(u"tm_in_nr_space");
    
    _TM_forward_sptr = std::make_shared<AffineTransformation<dataType> >();
    for (unsigned i=0; i<4; ++i)
        for (unsigned j=0; j<4; ++j)
            (*_TM_forward_sptr)[i][j] = float(tm[i][j]);
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
