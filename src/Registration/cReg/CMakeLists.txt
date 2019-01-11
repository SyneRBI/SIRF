#========================================================================
# Author: Richard Brown
# Copyright 2018 - 2019 University College London
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

# Add in the common files
SET(SOURCES "Registration.cpp" "NiftyRegistration.cpp" "NiftyAladinSym.cpp" "NiftyF3dSym.cpp" "ImageWeightedMean.cpp" "Resample.cpp" "NiftyResample.cpp" "Transformation.cpp" "AffineTransformation.cpp")

# Add in nifti_image wrappers
SET(SOURCES ${SOURCES} "NiftiImageData.cpp" "NiftiImageData3D.cpp" "NiftiImageData3DTensor.cpp" "NiftiImageData3DDeformation.cpp" "NiftiImageData3DDisplacement.cpp")

# If we're also wrapping to python or matlab, include the c-files
IF(BUILD_PYTHON OR BUILD_MATLAB)
	SET(SOURCES ${SOURCES} "cReg.cpp" "cReg_p.cpp")
ENDIF()

# Create library
ADD_LIBRARY(Reg ${SOURCES})

# Link to SIRF stuff
TARGET_LINK_LIBRARIES(Reg PUBLIC csirf)
target_link_libraries(Reg PUBLIC iutilities)
TARGET_INCLUDE_DIRECTORIES(Reg PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>$<INSTALL_INTERFACE:include/>"
)

find_package(OpenMP)

# NiftyReg
TARGET_LINK_LIBRARIES(Reg PUBLIC ${NiftyReg_Libs})
if (OpenMP_CXX_FOUND)
  TARGET_LINK_LIBRARIES(Reg PUBLIC OpenMP::OpenMP_CXX)
endif()

# Add libz (should be in /usr/local)
FIND_LIBRARY(LIBZ  z)
TARGET_LINK_LIBRARIES(Reg PUBLIC ${LIBZ})

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

# Install
INSTALL(TARGETS Reg  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib    )

# Need to include the directory to be able to compile the example programs
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

SET(REG_executables
    sirf_affine_to_disp
    sirf_print_nifti_info
    sirf_tensor_split_join
    sirf_change_datatype
    sirf_crop_image
    )

FOREACH(elem ${REG_executables})
    ADD_EXECUTABLE(${elem} ${elem}.cpp)
    TARGET_LINK_LIBRARIES(${elem} LINK_PUBLIC Reg)
    INSTALL(TARGETS ${elem} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
ENDFOREACH(elem ${REG_executables})

ADD_SUBDIRECTORY(tests)
