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

SET(NiftyReg_requiredLibs
    _reg_aladin
    _reg_blockMatching
    _reg_f3d
    _reg_maths
    _reg_ReadWriteImage
    _reg_resampling
    _reg_tools
    reg_nifti
    reg_png
    )
IF (NR_VERSION LESS 1.4)
    SET(NiftyReg_requiredLibs
        ${NiftyReg_requiredLibs}
        _reg_KLdivergence
        _reg_ReadWriteImage
        _reg_femTransformation
        _reg_globalTransformation
        _reg_localTransformation
        _reg_mutualinformation
        _reg_ssd
        _reg_thinPlateSpline
        reg_nrrd
    )
ELSEIF(NR_VERSION LESS 1.6)
    SET(NiftyReg_requiredLibs
        ${NiftyReg_requiredLibs}
        _reg_femTrans
        _reg_globalTrans
        _reg_localTrans
        _reg_measure
        )
ENDIF()
 
# Loop over each library, find it and add it. Then link
FOREACH(elem ${NiftyReg_requiredLibs})
    FIND_LIBRARY(NiftyReg_${elem} ${elem} HINTS ${NiftyReg_Install_DIR}/lib)
    SET(NiftyReg_Libs ${NiftyReg_${elem}} ${NiftyReg_Libs})
ENDFOREACH(elem ${NiftyReg_requiredLibs})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${NiftyReg_Libs})

# Add libz (should be in /usr/local)
FIND_LIBRARY(LIBZ  z)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBZ})

# Install
INSTALL(FILES   ${HEADERS}      DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib    )

# Need to include the directory to be able to compile the example programs
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(affine_to_disp affine_to_disp.cpp)
TARGET_LINK_LIBRARIES(affine_to_disp LINK_PUBLIC ${PROJECT_NAME})
INSTALL(TARGETS affine_to_disp DESTINATION ${CMAKE_INSTALL_PREFIX}/bin    )

ADD_SUBDIRECTORY(tests)
