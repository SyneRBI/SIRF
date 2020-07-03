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
\brief Base class for deformation models
\author Richard Brown
\author Kris Thielemans
\author SyneRBI
*/

#pragma once

namespace sirf {

// Forward declarations
template<class dataType> class NiftiImageData;
template<class dataType> class NiftiImageData3D;
template<class dataType> class NiftiImageData3DDeformation;

/*!
\ingroup Registration
\brief Base class for deformation models

\author Richard Brown
\author Kris Thielemans
\author SyneRBI
*/
 
template <class dataType>
   DeformationModel
  {
    virtual set_up(const NiftiImageData<dataType> &ref)
    {
      _ref = ref;
    }

    virtual void set_input(const Transformationparameters&pars)
    {
#if 0
      // if we add `Transformation::set()`
      this->_transformation_sptr->set(pars);
#else      
      this->_transformation_sptr = make_shared<transformType>(pars);
#endif      
    }
    virtual
      NiftiImageData3DDeformation<dataType>
      forward() const
    {
      return _transformation_sptr->get_as_deformation_field(ref);
    }

    std::unique<TransformationParameters>
      backward(const NiftiImageData3DDeformation<dataType>&)
      {
        return _transformation_sptr.back_from_deformation_field(ref);
      }

    virtual
      NiftiImageData3DDeformation<dataType>
      forward(const TransformationParameters& pars) const
    {
      set_input(pars);
      return forward();
    }

    std::unique_ptr<TransformationParameters>
      backward(const NiftiImageData3DDeformation<dataType>& dvf, const TransformationParameters&pars )
      {
        set_input(pars);
        return backward(dvf);
      }
  private:
    shared_ptr<Transformation> _transformation_sptr;
    NiftiImageData<dataType> _ref;
  };

}
