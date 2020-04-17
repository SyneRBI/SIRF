#========================================================================
# Author: Kris Thielemans
# Copyright 2016 - 2018 University College London
# Copyright 2017, 2018 Science Technology Facilities Council
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

set(CMAKE_POSITION_INDEPENDENT_CODE True)

# Check for existence of HKEM code
# TODO: replace if (EXISTS) with if (STIR version >= 3.2)
if(EXISTS "${STIR_INCLUDE_DIRS}/stir/KOSMAPOSL/KOSMAPOSLReconstruction.h")
  add_definitions(-DUSE_HKEM)
else()
  remove_definitions(-DUSE_HKEM)
endif()

# Check for existence of ListModeData.h
if (EXISTS "${STIR_INCLUDE_DIRS}/stir/listmode/ListModeData.h")
  add_definitions(-DSTIR_USE_LISTMODEDATA)
else()
  remove_definitions(-DSTIR_USE_LISTMODEDATA)
endif()

add_library(cstir 
    cstir_p.cpp cstir_tw.cpp stir_data_containers.cpp stir_x.cpp cstir.cpp)
target_include_directories(cstir PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>$<INSTALL_INTERFACE:include>")
target_include_directories(cstir PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")
target_include_directories(cstir PUBLIC "${STIR_INCLUDE_DIRS}")

target_link_libraries(cstir csirf iutilities)
target_link_libraries(cstir "${STIR_LIBRARIES}")
# Add boost library dependencies
if((CMAKE_VERSION VERSION_LESS 3.5.0) OR (NOT _Boost_IMPORTED_TARGETS))
  # This is harder than it should be on older CMake versions to be able to cope with
  # spaces in filenames.
  foreach(C SYSTEM FILESYSTEM THREAD DATE_TIME CHRONO)
    target_link_libraries(cstir optimized "${Boost_${C}_LIBRARY_RELEASE}")
    target_link_libraries(cstir debug  "${Boost_${C}_LIBRARY_DEBUG}")
  endforeach()
else()
  # Nice and simple for recent CMake (which knows about your Boost version)
  target_link_libraries(cstir Boost::system Boost::filesystem Boost::thread Boost::date_time Boost::chrono)
endif()

ADD_SUBDIRECTORY(tests)
