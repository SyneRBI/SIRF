/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

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

#include "sirf/Reg/SPMRegistration.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include <sys/stat.h>
#include <MatlabEngine.hpp>
#include <boost/filesystem.hpp>
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include <codecvt>

using namespace sirf;
using namespace matlab::engine;
using namespace matlab::data;

inline bool check_file_exists(const std::string& filename, const bool existance_allowed) {
    struct stat buffer;
    const bool file_exists (stat (filename.c_str(), &buffer) == 0);
    if (file_exists && !existance_allowed)
        throw std::runtime_error("SPMRegistration<dataType>::process(): file exists: " + filename);
    return file_exists;
}

template<class dataType>
void SPMRegistration<dataType>::set_working_folder(const std::string &working_folder)
{
    // Make sure it's absolute
    _working_folder = boost::filesystem::absolute(working_folder).string();
}

template<class dataType>
void SPMRegistration<dataType>::process()
{
    // Check the paramters that are NOT set via the parameter file have been set.
    this->check_parameters();

    // If the working folder doesn't already exist, and delete is desired, remember to delete it
    bool need_to_delete_working_folder = false;
    if (_delete_temp_files && !check_file_exists(_working_folder,true))
        need_to_delete_working_folder = true;

    const size_t num_flo_ims = this->_floating_images.size();

    // Reference
    const std::string ref_filename = _working_folder + "/ref.nii";
    check_file_exists(ref_filename, _working_folder_overwrite);
    NiftiBasedRegistration<dataType>::convert_to_NiftiImageData_if_not_already(this->_reference_image_nifti_sptr, this->_reference_image_sptr);
    this->_reference_image_nifti_sptr->write(ref_filename);

    // Floatings
    std::vector<std::string> flo_filenames(num_flo_ims), resliced_filenames(num_flo_ims);
    this->_floating_images_nifti.resize(num_flo_ims);
    for (unsigned i=0; i<num_flo_ims; ++i) {
        flo_filenames.at(i) = _working_folder + "/flo" + std::to_string(i) + ".nii";
        resliced_filenames.at(i) = _working_folder + "/rflo" + std::to_string(i) + ".nii";
        check_file_exists(flo_filenames.at(i), _working_folder_overwrite);
        check_file_exists(resliced_filenames.at(i), _working_folder_overwrite);
        NiftiBasedRegistration<dataType>::convert_to_NiftiImageData_if_not_already(this->_floating_images_nifti.at(i), this->_floating_images.at(i));
        this->_floating_images_nifti.at(i)->write(flo_filenames.at(i));
    }

    // Start MATLAB engine synchronously
    if (!_matlab_uptr) {
        _matlab_uptr = startMATLAB(std::vector<String>({u"-nojvm"}));
        std::cout << "Started MATLAB Engine" << std::endl;
    }

    // We'll need to be able to convert commands from std::string to std::u16string
    std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t> convertor;

    // Call spm_realign
    std::string spm_realign_command = "spm_output = spm_realign(char('" + ref_filename;
    for (unsigned i=0; i<num_flo_ims; ++i)
        spm_realign_command +=  + "','" + flo_filenames.at(i);
    spm_realign_command += "'),struct('quality',1));";
    _matlab_uptr->eval(convertor.from_bytes(spm_realign_command));

    // Call spm_reslice
    const std::string spm_reslice_command = "spm_reslice(spm_output,struct('which',[1,0]));";
    _matlab_uptr->eval(convertor.from_bytes(spm_reslice_command));

    // Read warped image back in
    this->_warped_images_nifti.resize(num_flo_ims);
    for (unsigned i=0; i<num_flo_ims; ++i)
        this->_warped_images_nifti.at(i) = std::make_shared<NiftiImageData3D<dataType> >(resliced_filenames.at(i));

    // Clean up
    if (_delete_temp_files) {
        remove( ref_filename.c_str() );
        for (unsigned i=0; i<num_flo_ims; ++i) {
            remove( flo_filenames.at(i).c_str() );
            remove( resliced_filenames.at(i).c_str() );
        }
        if (need_to_delete_working_folder)
            remove( _working_folder.c_str() );
    }

    // Get the transformation matrices
    _TMs_fwd.resize(num_flo_ims);
    _TMs_inv.resize(num_flo_ims);
    for (unsigned a=0; a<num_flo_ims; ++a) {
        const std::string tm_in_nr_space_command = "tm_in_nr_space = inv(spm_output(" + std::to_string(a+2) + ").mat / spm_output(" + std::to_string(a+2) + ").private.mat0);";
        _matlab_uptr->eval(convertor.from_bytes(tm_in_nr_space_command));
        TypedArray<double> tm = _matlab_uptr->getVariable(u"tm_in_nr_space");

        _TMs_fwd.at(a) = std::make_shared<AffineTransformation<dataType> >();
        for (unsigned i=0; i<4; ++i)
            for (unsigned j=0; j<4; ++j)
                (*_TMs_fwd.at(a))[i][j] = float(tm[i][j]);
        _TMs_inv.at(a) = std::make_shared<AffineTransformation<dataType> >(_TMs_fwd.at(a)->get_inverse());
    }

    // Get as deformation and displacement
    this->_disp_fwd_images.resize(num_flo_ims);
    this->_disp_inv_images.resize(num_flo_ims);
    for (unsigned i=0; i<num_flo_ims; ++i) {
        NiftiImageData3DDeformation<dataType> def_fwd = _TMs_fwd.at(i)->get_as_deformation_field(*this->_reference_image_nifti_sptr);
        NiftiImageData3DDeformation<dataType> def_inv = *def_fwd.get_inverse(this->_floating_images_nifti.at(i));
        this->_disp_fwd_images.at(i) = std::make_shared<NiftiImageData3DDisplacement<dataType> >(def_fwd);
        this->_disp_inv_images.at(i) = std::make_shared<NiftiImageData3DDisplacement<dataType> >(def_inv);
    }
    // The output should be a clone of the reference image, with data filled in from the nifti image
    this->_warped_images.resize(num_flo_ims);
    for (unsigned i=0; i<num_flo_ims; ++i) {
        this->_warped_images.at(i) = this->_reference_image_sptr->clone();
        this->_warped_images.at(i)->fill(*this->_warped_images_nifti.at(i));
    }
}

template<class dataType>
void SPMRegistration<dataType>::check_parameters() const
{
    // Call base class
    Registration<dataType>::check_parameters();

    if (_working_folder.empty())
        throw std::runtime_error("SPMRegistration<dataType>::check_parameters(): Missing working folder.");
}

namespace sirf {
template class SPMRegistration<float>;
}
