#========================================================================
# Author: Kris Thielemans
# Copyright 2016, 2017, 2019 University College London
# Copyright 2018 Science Technology Facilities Council
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
# Commented out, as this is now done at the top-level
#find_package(Boost COMPONENTS system filesystem thread date_time chrono REQUIRED)
# For Visual Studio we have to disable the auto-linking feature of boost
# where just including a boost file automatically adds it to the linker path.
# Although this sounds great, it sadly breaks because of conflicting file-paths when linking etc etc.
# In any case, we need to add the libraries by hand for other systems.
# See http://stackoverflow.com/questions/32252016/cmake-visual-studio-build-looks-for-wrong-library
#add_definitions(-DBOOST_ALL_NO_LIB)

if (SIRF_INSTALL_DEPENDENCIES AND WIN32)
    set(Boost_DLL_DIR ${Boost_LIBRARY_DIR_RELEASE})
    message(STATUS "Install boost shared libraries from ${Boost_DLL_DIR} ")
	foreach (__boost_lib system filesystem thread date_time chrono program_options)
	   file(GLOB Boost_DLL "${Boost_DLL_DIR}/boost_${__boost_lib}*.dll")
       install( FILES ${Boost_DLL} DESTINATION bin )
	endforeach()
  endif()
	
add_library(cgadgetron cgadgetron.cpp gadgetron_x.cpp gadgetron_data_containers.cpp gadgetron_client.cpp 
            gadgetron_fftw.cpp ismrmrd_fftw.cpp ismrmrd_phantom.cpp shepp_logan_phantom.cpp)

set (cGadgetron_INCLUDE_DIR "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>$<INSTALL_INTERFACE:include>")
target_include_directories(cgadgetron PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")
# copy to parent scope
#set (cGadgetron_INCLUDE_DIR "${cGadgetron_INCLUDE_DIR}" PARENT_SCOPE)

target_include_directories(cgadgetron PUBLIC "${cGadgetron_INCLUDE_DIR}")
target_include_directories(cgadgetron PRIVATE "${FFTW3_INCLUDE_DIR}")
target_include_directories(cgadgetron PUBLIC "${ISMRMRD_INCLUDE_DIR}")

target_link_libraries(cgadgetron iutilities csirf)
# Add boost library dependencies
if((CMAKE_VERSION VERSION_LESS 3.5.0) OR (NOT _Boost_IMPORTED_TARGETS))
  # This is harder than it should be on older CMake versions to be able to cope with
  # spaces in filenames.
  foreach(C SYSTEM FILESYSTEM THREAD DATE_TIME CHRONO)
    target_link_libraries(cgadgetron optimized "${Boost_${C}_LIBRARY_RELEASE}")
    target_link_libraries(cgadgetron debug  "${Boost_${C}_LIBRARY_DEBUG}")
  endforeach()
else()
  # Nice and simple for recent CMake (which knows about your Boost version)
  target_link_libraries(cgadgetron Boost::system Boost::filesystem Boost::thread Boost::date_time Boost::chrono)
endif()

# Note: cannot use ISMRMRD_LIBRARIES on Windows as it generally contains 
# a list of filenames with spaces. There doesn't seem to be a way to pass this through (strange).
# Luckily, we know what libraries it uses
find_library(ismrmrd_full_path "ismrmrd" "${ISMRMRD_LIBRARY_DIRS}")
if(NOT ismrmrd_full_path)
  message(FATAL_ERROR "ismrmrd not found")
endif()
target_link_libraries(cgadgetron ${ismrmrd_full_path})
target_link_libraries(cgadgetron "${FFTW3_LIBRARIES}")
target_link_libraries(cgadgetron "${HDF5_LIBRARIES}")

ADD_SUBDIRECTORY(tests)
