#========================================================================
# Author: Richard Brown
# Copyright 2018 - 2020 University College London
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

set(CMAKE_POSITION_INDEPENDENT_CODE True)

# Configure the aladin/f3d files dependent on NIFTYREG version number
include(${CMAKE_CURRENT_SOURCE_DIR}/../callable_NR_methods.cmake)
get_parser_methods("../callable_NR_aladin_methods.txt" "reg_aladin_sym" "${NIFTYREG_VERSION}" NR_aladin_parser_methods)
get_parser_methods("../callable_NR_f3d_methods.txt"    "reg_f3d_sym"    "${NIFTYREG_VERSION}" NR_f3d_parser_methods   )
get_runtime_methods("../callable_NR_aladin_methods.txt" "${NIFTYREG_VERSION}" NR_aladin_runtime_methods)
get_runtime_methods("../callable_NR_f3d_methods.txt"    "${NIFTYREG_VERSION}" NR_f3d_runtime_methods   )
get_list_of_wrapped_methods("../callable_NR_aladin_methods.txt" "${NIFTYREG_VERSION}" NR_aladin_list_methods)
get_list_of_wrapped_methods("../callable_NR_f3d_methods.txt"    "${NIFTYREG_VERSION}" NR_f3d_list_methods   )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/NiftyAladinSym.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/NiftyAladinSym.cpp")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/NiftyF3dSym.cpp.in"    "${CMAKE_CURRENT_BINARY_DIR}/NiftyF3dSym.cpp"   )

# Add in the common files
SET(SOURCES
  "Registration.cpp"
  "NiftiBasedRegistration.cpp"
  "NiftyRegistration.cpp"
  "${CMAKE_CURRENT_BINARY_DIR}/NiftyAladinSym.cpp"
  "${CMAKE_CURRENT_BINARY_DIR}/NiftyF3dSym.cpp"
  "ImageWeightedMean.cpp"
  "Resample.cpp"
  "NiftyResample.cpp"
  "Transformation.cpp"
  "AffineTransformation.cpp"
  "Quaternion.cpp"
  "NiftiImageData.cpp"
  "NiftiImageData3D.cpp"
  "NiftiImageData3DTensor.cpp"
  "NiftiImageData3DDeformation.cpp"
  "NiftiImageData3DDisplacement.cpp")

# If we're also wrapping to python or matlab, include the c-files
IF(BUILD_PYTHON OR BUILD_MATLAB)
	SET(SOURCES ${SOURCES} "cReg.cpp" "cReg_p.cpp")
ENDIF()

# Add SPM source
IF(SPM_FOUND)
  SET(SOURCES ${SOURCES} "SPMRegistration.cpp")
ENDIF()

# Create library
ADD_LIBRARY(Reg ${SOURCES})

# Link to SIRF stuff
TARGET_LINK_LIBRARIES(Reg PUBLIC csirf)
target_link_libraries(Reg PUBLIC iutilities)
TARGET_INCLUDE_DIRECTORIES(Reg PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>$<INSTALL_INTERFACE:include/>"
)

# NIFTYREG
TARGET_LINK_LIBRARIES(Reg PUBLIC ${NR_libs_full_path} ${NR_libs_full_path})
TARGET_INCLUDE_DIRECTORIES(Reg BEFORE PUBLIC "${NIFTYREG_INCLUDE_DIRS}")

# link against NiftyMoMo
TARGET_LINK_LIBRARIES(Reg PUBLIC NiftyMoMo)

# Temporary fix for linking CUDA libraries if necessary
if (NIFTYREG_BUILT_WITH_CUDA)
  find_package(CUDA REQUIRED)
  target_link_libraries(Reg PUBLIC ${CUDA_CUDART_LIBRARY})
endif()

# Add boost library dependencies
if((CMAKE_VERSION VERSION_LESS 3.5.0) OR (NOT _Boost_IMPORTED_TARGETS))
  # This is harder than it should be on older CMake versions to be able to cope with
  # spaces in filenames.
  foreach(C SYSTEM FILESYSTEM)
    target_link_libraries(Reg PUBLIC optimized "${Boost_${C}_LIBRARY_RELEASE}")
    target_link_libraries(Reg PUBLIC debug  "${Boost_${C}_LIBRARY_DEBUG}")
  endforeach()
else()
  # Nice and simple for recent CMake (which knows about your Boost version)
  target_link_libraries(Reg PUBLIC Boost::system Boost::filesystem)
endif()

# Link Matlab for SPM
IF(SPM_FOUND)
  TARGET_INCLUDE_DIRECTORIES(Reg PUBLIC "${Matlab_INCLUDE_DIRS}")
  TARGET_LINK_LIBRARIES(Reg PUBLIC "${Matlab_LIBRARIES}")
  target_compile_definitions(Reg PUBLIC SIRF_SPM)
ENDIF()

# Install
INSTALL(TARGETS Reg  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

SET(REG_executables
    sirf_affine_to_disp
    sirf_print_nifti_info
    sirf_tensor_split_join
    sirf_change_datatype
    sirf_crop_image
    sirf_flip_or_mirror_image
    sirf_nifti_maths
    sirf_inverse_transform
    )

FOREACH(elem ${REG_executables})
    ADD_EXECUTABLE(${elem} ${elem}.cpp)
    TARGET_LINK_LIBRARIES(${elem} LINK_PUBLIC Reg)
    INSTALL(TARGETS ${elem} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
ENDFOREACH(elem ${REG_executables})

ADD_SUBDIRECTORY(tests)
