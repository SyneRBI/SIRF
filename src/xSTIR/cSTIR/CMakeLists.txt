#========================================================================
# Author: Kris Thielemans
# Author: Richard Brown
# Author: Dimitra Kyriakopoulou
# Copyright 2016 - 2020 University College London
# Copyright 2026 University College London
# Copyright 2026 Biomedical Research Foundation, Academy of Athens
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

include(CheckIncludeFileCXX)

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

if (TARGET CUDA::cudart)
  get_target_property(_cudart_include_dirs CUDA::cudart INTERFACE_INCLUDE_DIRECTORIES)
  set(_saved_required_includes "${CMAKE_REQUIRED_INCLUDES}")
  if (_cudart_include_dirs)
    set(CMAKE_REQUIRED_INCLUDES ${_cudart_include_dirs})
  endif()
  check_include_file_cxx(cuda_runtime_api.h HAS_CUDA_RUNTIME_API)
  set(CMAKE_REQUIRED_INCLUDES "${_saved_required_includes}")
  if (HAS_CUDA_RUNTIME_API)
    target_compile_definitions(cstir PUBLIC HAS_CUDA_RUNTIME_API)
    target_link_libraries(cstir CUDA::cudart)
  endif()
endif()

target_link_libraries(cstir csirf iutilities)
target_link_libraries(cstir "${STIR_LIBRARIES}")
# Add boost library dependencies
target_link_libraries(cstir Boost::system Boost::filesystem Boost::thread Boost::date_time Boost::chrono)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  COMPONENT Development
  DESTINATION include)
install(TARGETS cstir EXPORT SIRFTargets
  COMPONENT Development
  DESTINATION lib)

ADD_SUBDIRECTORY(tests)
