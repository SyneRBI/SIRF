/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

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
\ingroup Registration
\brief Wrapper around SPM's registration class.

\author Richard Brown
\author SyneRBI
*/

#pragma once

#include "sirf/Reg/NiftiBasedRegistration.h"
#include <MatlabEngine.hpp>

namespace sirf {

/// Forward declarations
template<class dataType> class AffineTransformation;

/*!
\ingroup Registration
\brief Wrapper around SPM's registration class.

\author Richard Brown
\author SyneRBI
*/
template<class dataType> class SPMRegistration : public NiftiBasedRegistration<dataType>
{
public:

    /// Destructor
    virtual ~SPMRegistration();

    /// Process
    void process();

    /// Set working folder
    void set_working_folder(const std::string &working_folder);

    /// Set file overwrite in working folder
    void set_working_folder_file_overwrite(const bool working_folder_overwrite) { _working_folder_overwrite = working_folder_overwrite; }

    /// Delete temporary files
    void set_delete_temp_files(const bool delete_temp_files) { _delete_temp_files = delete_temp_files; }

    /// Get registered image
    virtual const std::shared_ptr<const ImageData> get_output_sptr(const unsigned idx = 0) const;

    /// Get forwards transformation matrix
    virtual const std::shared_ptr<const AffineTransformation<float> > get_transformation_matrix_forward_sptr(const unsigned idx = 0) const { return _TMs_fwd.at(idx); }

    /// Get inverse transformation matrix
    virtual const std::shared_ptr<const AffineTransformation<float> > get_transformation_matrix_inverse_sptr(const unsigned idx = 0) const { return _TMs_inv.at(idx); }

    /// Get forward deformation field image
    virtual const std::shared_ptr<const Transformation<dataType> > get_deformation_field_forward_sptr(const unsigned idx = 0) const;

    /// Get inverse deformation field image
    virtual const std::shared_ptr<const Transformation<dataType> > get_deformation_field_inverse_sptr(const unsigned idx = 0) const;

protected:

    /// Parse parameter file
    virtual void parse_parameter_file() {}

    /// Check parameters
    virtual void check_parameters() const;

    /// Set any extra parameters
    virtual void set_parameters() {}

    /// working folder
    std::string _working_folder = "";
    /// Overwrite files already in working folder
    bool _working_folder_overwrite = false;
    /// Delete temp files
    bool _delete_temp_files = false;

    /// Forwards transformation matrix
    std::vector<std::shared_ptr<AffineTransformation<float> > > _TMs_fwd;
    /// Inverse transformation matrix
    std::vector<std::shared_ptr<AffineTransformation<float> > > _TMs_inv;
    /// Matlab instance
    std::unique_ptr<matlab::engine::MATLABEngine> _matlab_uptr;

    /// Filenames to potentially delete
    std::vector<std::string> _filenames_to_delete;
    /// Folders to potentially delete
    std::vector<std::string> _folders_to_delete;
    /// Resliced images
    std::vector<std::string> _resliced_filenames;

};
}
