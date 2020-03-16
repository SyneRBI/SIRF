
#========================================================================
# Author: Richard Brown, Kris Thielemans
# Copyright 2016 - 2019 University College London
# Copyright 2016 - 2018 Science Technology Facilities Council
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0.txt
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#=========================================================================

ADD_SUBDIRECTORY(iUtilities)
ADD_SUBDIRECTORY(common)


##########################################################################
#                              Gadgetron                                 #
##########################################################################
option(DISABLE_Gadgetron "Disable building the SIRF interface to Gadgetron" OFF)
if (DISABLE_Gadgetron)
  message(STATUS "Gadgetron support disabled.")
else()
  # include explicit search for HDF5 for older ISMRMRD versions.
  find_package(HDF5 1.8 COMPONENTS C REQUIRED)
  find_package(ISMRMRD REQUIRED)
  # Add ISMRMRD to search path for FFTW3
  set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${ISMRMRD_DIR}")
  find_package(FFTW3 COMPONENTS single REQUIRED)
  ADD_SUBDIRECTORY(xGadgetron)
endif()



##########################################################################
#                             Registration                               #
##########################################################################
option(DISABLE_Registration "Disable building the SIRF registration package" OFF)
if (DISABLE_Registration)
  message(STATUS "Registration support disabled.")
else()
  FIND_PACKAGE(NIFTYREG 1.5.61 REQUIRED)
  # NIFTYREG
  FOREACH(NR_lib ${NIFTYREG_LIBRARIES})
    find_library(${NR_lib}_full_path "${NR_lib}" "${NIFTYREG_LIBRARY_DIRS}")
    if(NOT ${NR_lib}_full_path)
      message(FATAL_ERROR "${NR_lib} not found")
    endif()
    SET(NR_libs_full_path "${NR_libs_full_path};${${NR_lib}_full_path}")
  ENDFOREACH()
  # If niftyreg was bulit with OpenMP
  if (NIFTYREG_BUILT_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    if (OpenMP_CXX_FOUND)
      SET(NR_libs_full_path "${NR_libs_full_path};OpenMP::OpenMP_CXX")
    endif()
  endif()
  # If niftyreg was bulit with CUDA
  if (NIFTYREG_BUILT_WITH_CUDA)
    find_package(CUDA REQUIRED)
    SET(NR_libs_full_path "${NR_libs_full_path};${CUDA_CUDA_LIBRARY};${CUDA_CUDART_LIBRARY}")
  endif()
  ADD_SUBDIRECTORY(Registration)
endif()



##########################################################################
#                                 STIR                                   #
##########################################################################
option(DISABLE_STIR "Disable building the SIRF interface to STIR" OFF)
if (DISABLE_STIR)
  message(STATUS "STIR support disabled.")
else()
  find_package(STIR 4.0.0 REQUIRED)
  if (STIR_WITH_NIFTYPET_PROJECTOR)
    set(NIFTYPET_BOOL_STR "1")
  	message(STATUS "STIR was built with NiftyPET, GPU projectors will be enabled.")
    add_compile_definitions(STIR_WITH_NIFTYPET_PROJECTOR)
  else()
    MESSAGE(STATUS "STIR not built with NiftyPET.")
    set(NIFTYPET_BOOL_STR "0")
  endif()
  ADD_SUBDIRECTORY(xSTIR)
endif()


##########################################################################
#                             Synergistic                                #
##########################################################################
# Synergistic code currently only works if STIR (with ITK), ISMRMRD and Registration are built
option(DISABLE_Synergistic "Disable building the synergistic code" OFF)
if (DISABLE_Synergistic)
  message(STATUS "Synergistic code disabled.")
else()
  if ((NOT DISABLE_STIR) AND (NOT DISABLE_Gadgetron) AND (NOT DISABLE_Registration) AND "${STIR_BUILT_WITH_ITK}")
    MESSAGE(STATUS "Registration, ISMRMRD and STIR (with ITK) have been built. Building synergistic code.")
    ADD_SUBDIRECTORY(Synergistic)
  else()
    MESSAGE(STATUS "One or more of the following components are missing: Registration, ISMRMRD and STIR (with ITK). Synergistic code will not be built.")
    if (DISABLE_STIR)
      MESSAGE(STATUS "STIR missing.")
    elseif(NOT STIR_BUILT_WITH_ITK)
      MESSAGE(STATUS "STIR not built with ITK.")
    endif()
    if (DISABLE_Gadgetron)
      MESSAGE(STATUS "ISMRMRD missing.")
    endif()
    if (DISABLE_Registration)
      MESSAGE(STATUS "Registration missing.")
    endif()
  endif()
endif()