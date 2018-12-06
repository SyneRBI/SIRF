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

set(CMAKE_POSITION_INDEPENDENT_CODE True)

# Add in the common files
SET(SOURCES "SIRFReg.cpp" "SIRFRegNiftyRegistration.cpp" "SIRFRegNiftyAladinSym.cpp" "SIRFRegNiftyF3dSym.cpp" "SIRFRegImageWeightedMean.cpp" "SIRFRegResample.cpp" "SIRFRegNiftyResample.cpp" "SIRFRegTransformation.cpp" "SIRFRegAffineTransformation.cpp")
SET(HEADERS  "SIRFReg.h"    "SIRFRegNiftyRegistration.h"  "SIRFRegNiftyAladinSym.h"   "SIRFRegNiftyF3dSym.h"   "SIRFRegImageWeightedMean.h"   "SIRFRegResample.h"   "SIRFRegNiftyResample.h"   "SIRFRegTransformation.h"   "SIRFRegAffineTransformation.h" "SIRFRegParser.h"   "SIRFRegParserKey.h" )

# Add in nifti_image wrappers
SET(SOURCES ${SOURCES} "NiftiImageData.cpp" "NiftiImageData3D.cpp" "NiftiImageData3DTensor.cpp" "NiftiImageData3DDeformation.cpp" "NiftiImageData3DDisplacement.cpp")
SET(HEADERS ${HEADERS}  "NiftiImageData.h"   "NiftiImageData3D.h"   "NiftiImageData3DTensor.h"   "NiftiImageData3DDeformation.h"   "NiftiImageData3DDisplacement.h")

# If we're also wrapping to python or matlab, include the c-files
IF(BUILD_PYTHON OR BUILD_MATLAB)
	SET(SOURCES ${SOURCES} "csirfreg.cpp" "csirfreg_p.cpp")
	SET(HEADERS ${HEADERS}  "csirfreg.h"   "csirfreg_p.h" )
ENDIF()

# Create library
ADD_LIBRARY(SIRFReg ${SOURCES})

# Link to SIRF stuff
TARGET_LINK_LIBRARIES(SIRFReg common)
target_link_libraries(SIRFReg iutilities)

# NiftyReg
TARGET_LINK_LIBRARIES(SIRFReg ${NiftyReg_Libs})

# Add libz (should be in /usr/local)
FIND_LIBRARY(LIBZ  z)
TARGET_LINK_LIBRARIES(SIRFReg ${LIBZ})

# Add boost library dependencies
if((CMAKE_VERSION VERSION_LESS 3.5.0) OR (NOT _Boost_IMPORTED_TARGETS))
  # This is harder than it should be on older CMake versions to be able to cope with
  # spaces in filenames.
  foreach(C SYSTEM FILESYSTEM)
    target_link_libraries(SIRFReg optimized "${Boost_${C}_LIBRARY_RELEASE}")
    target_link_libraries(SIRFReg debug  "${Boost_${C}_LIBRARY_DEBUG}")
  endforeach()
else()
  # Nice and simple for recent CMake (which knows about your Boost version)
  target_link_libraries(SIRFReg Boost::system Boost::filesystem)
endif()

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
