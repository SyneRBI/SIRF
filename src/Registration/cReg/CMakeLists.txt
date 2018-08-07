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

# Set the dependencies list
SET(${proj}_DEPENDENCIES "nifti;Boost")

# Add in the common files
SET(SOURCES "SIRFReg.cpp" "SIRFRegNiftyAladinSym.cpp" "SIRFRegNiftyF3dSym.cpp" "SIRFRegImageWeightedMean.cpp" "SIRFRegNiftyResample.cpp" "SIRFRegMisc.cpp" "SIRFRegActivityCorrect.cpp" "SIRFRegExtractTAC.cpp" "SIRFImageData.cpp")
SET(HEADERS  "SIRFReg.h"   "SIRFRegNiftyAladinSym.h"   "SIRFRegNiftyF3dSym.h"   "SIRFRegImageWeightedMean.h"   "SIRFRegNiftyResample.h"   "SIRFRegMisc.h"   "SIRFRegActivityCorrect.h"   "SIRFRegExtractTAC.h"   "SIRFImageData.h" "SIRFRegParser.h"   "SIRFRegParserKey.h" )

# Create library
ADD_LIBRARY(${PROJECT_NAME} SHARED ${SOURCES})

# Compile against STIR
IF(STIR_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} cstir)
ENDIF(STIR_FOUND)

# Boost
FIND_PACKAGE(Boost 1.63.0 COMPONENTS filesystem)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES})

# NiftyReg
TARGET_LINK_LIBRARIES(SIRFReg ${NiftyReg_Libs})

# Add libz (should be in /usr/local)
FIND_LIBRARY(LIBZ  z)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBZ})

# Install
INSTALL(FILES   ${HEADERS}      DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib    )

# Need to include the directory to be able to compile the example programs
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SIRFREG_executables
    affine_to_disp
    sirfreg_aladin)

FOREACH(elem ${SIRFREG_executables})
    ADD_EXECUTABLE(${elem} ${elem}.cpp)
    TARGET_LINK_LIBRARIES(${elem} LINK_PUBLIC ${PROJECT_NAME})
    INSTALL(TARGETS ${elem} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
ENDFOREACH(elem ${SIRFREG_executables})

ADD_SUBDIRECTORY(tests)
