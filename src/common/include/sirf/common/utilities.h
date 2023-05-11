/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2023 Rutherford Appleton Laboratory STFC
Copyright 2023 University College London

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
\ingroup Common

\author Evgueni Ovtchinnikov
\author Kris Thielemans
*/
#ifndef SIRF_UTILITIES
#define SIRF_UTILITIES

#include <string>

namespace sirf {
    //! return the path-separator used by the OS
    /*!
      \ingroup Common
      Usually this will return a forward-slash, but it could be a backslash on Windows.
     */
    char path_separator();
    ///@{
    //! concatenate path strings
    /*!
      \ingroup Common
      \par Example:
      \code
      std::string filename = "somefile.txt";
      std:string full_path = append_path("/usr/local", "share", filename);
      \endcode
      \par warning
      The code does not check if there are already trailing path-separators,
      nor if the resulting path exists.
    */
    template <typename T>
      std::string append_path(std::string path, T a)
      {
        return path + path_separator() + std::string(a);
      }
    template <typename T, typename... Ts>
      std::string append_path(std::string path, T a, Ts... args)
    {
      return append_path(append_path(path, a), args...);
    }
    ///@{
    //! Returns where the examples are installed
    /*!
      \ingroup Common
      \par Example:
      \code
      std::string p = examples_data_path("PET");
      \endcode
    */
    std::string examples_data_path(const char* data_type);
}

#endif
