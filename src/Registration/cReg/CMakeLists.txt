#========================================================================
# Author: Richard Brown
# Copyright 2016, 2017 University College London
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

# Add in the common files
SET(SOURCES "SIRFReg.cpp" "SIRFRegNiftyAladinSym.cpp" "SIRFRegNiftyF3dSym.cpp" "SIRFRegImageWeightedMean.cpp" "SIRFRegNiftyResample.cpp" "SIRFRegMisc.cpp" "SIRFRegTransformation.cpp" "SIRFRegMat44.cpp")
SET(HEADERS  "SIRFReg.h"   "SIRFRegNiftyAladinSym.h"   "SIRFRegNiftyF3dSym.h"   "SIRFRegImageWeightedMean.h"  "SIRFRegNiftyResample.h"   "SIRFRegMisc.h"    "SIRFRegTransformation.h"   "SIRFRegMat44.h" "SIRFRegParser.h"   "SIRFRegParserKey.h" )

# Add in nifti_image wrappers
SET(SOURCES ${SOURCES} "NiftiImage.cpp" "NiftiImage3D.cpp" "NiftiImage3DTensor.cpp" "NiftiImage3DDeformation.cpp" "NiftiImage3DDisplacement.cpp")
SET(HEADERS ${HEADERS}  "NiftiImage.h"   "NiftiImage3D.h"   "NiftiImage3DTensor.h"   "NiftiImage3DDeformation.h"   "NiftiImage3DDisplacement.h")

# If we're also wrapping to python or matlab, include the c-files
IF(SIRFREG_BUILD_PYTHON OR SIRFREG_BUILD_MATLAB)
	SET(SOURCES ${SOURCES} "csirfreg.cpp" "csirfreg_p.cpp")
	SET(HEADERS ${HEADERS}  "csirfreg.h"   "csirfreg_p.h" )
ENDIF()

# Create library
ADD_LIBRARY(SIRFReg ${SOURCES})
#TARGET_LINK_LIBRARIES(SIRFReg SIRF_common)

# Link against stir and gadgetron
TARGET_LINK_LIBRARIES(SIRFReg cstir)
#TARGET_LINK_LIBRARIES(SIRFReg cgadgetron)

# Boost
FIND_PACKAGE(Boost 1.63.0 COMPONENTS filesystem)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(SIRFReg ${Boost_LIBRARIES})

# NiftyReg
TARGET_LINK_LIBRARIES(SIRFReg ${NiftyReg_Libs})

# Add libz (should be in /usr/local)
FIND_LIBRARY(LIBZ  z)
TARGET_LINK_LIBRARIES(SIRFReg ${LIBZ})

# Install
INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
INSTALL(TARGETS SIRFReg  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib    )

# Need to include the directory to be able to compile the example programs
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SIRFREG_executables
    sirfreg_affine_to_disp
    sirfreg_print_nifti_info
    sirfreg_aladin
    sirfreg_tensor_split_join
    sirfreg_change_datatype
    sirfreg_crop_image
    )

FOREACH(elem ${SIRFREG_executables})
    ADD_EXECUTABLE(${elem} ${elem}.cpp)
    TARGET_LINK_LIBRARIES(${elem} LINK_PUBLIC SIRFReg)
    INSTALL(TARGETS ${elem} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
ENDFOREACH(elem ${SIRFREG_executables})

ADD_SUBDIRECTORY(tests)
