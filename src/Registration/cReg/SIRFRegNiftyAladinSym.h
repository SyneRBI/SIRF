/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

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

#ifndef _SIRFREGNIFTYALADINSYM_H_
#define _SIRFREGNIFTYALADINSYM_H_


#include "SIRFReg.h"
#include "SIRFRegMat44.h"
#if NIFTYREG_VER_1_5
#include <_reg_aladin_sym.h>
#elif NIFTYREG_VER_1_3
#include <_reg_aladin.h>
#endif

namespace sirf {
/// Wrapper around NiftyReg's aladin class for rigid and affine transformations
template<class T> class SIRFRegNiftyAladinSym : public SIRFReg
{
public:

    /// Update
    void update();

    /// Get forwards transformation matrix
    const SIRFRegMat44 &get_transformation_matrix_fwrd() const { return _TM_fwrd; }

    /// Get backwards transformation matrix
    const SIRFRegMat44 &get_transformation_matrix_back() const { return _TM_back; }

protected:

    /// Parse parameter file
    virtual void parse_parameter_file();

    /// Register object
#if NIFTYREG_VER_1_5
    std::shared_ptr<reg_aladin_sym<T> > _registration_sptr;
#elif NIFTYREG_VER_1_3
    std::shared_ptr<reg_aladin<T> > _registration_sptr;
#endif

    /// Forwards transformation matrix
    SIRFRegMat44 _TM_fwrd;
    /// Backwards transformation matrix
    SIRFRegMat44 _TM_back;
};
}

#endif
