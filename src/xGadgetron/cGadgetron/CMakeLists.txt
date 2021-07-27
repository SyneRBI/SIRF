#========================================================================
# Author: Kris Thielemans
# Copyright 2016, 2017, 2019, 2020 University College London
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

set(CGADGETRON_SOURCES cgadgetron.cpp gadgetron_x.cpp gadgetron_data_containers.cpp gadgetron_client.cpp
    gadgetron_fftw.cpp TrajectoryPreparation.cpp FourierEncoding.cpp)
    

find_library( GT_CPUCORE NAMES gadgetron_toolbox_cpucore
                NAMES_PER_DIR
                PATHS ${CMAKE_INSTALL_PREFIX}
                )
find_library( GT_FFT NAMES gadgetron_toolbox_cpufft
                NAMES_PER_DIR
                PATHS ${CMAKE_INSTALL_PREFIX}
                )
find_library( GT_NFFT NAMES gadgetron_toolbox_cpunfft
                NAMES_PER_DIR
                PATHS ${CMAKE_INSTALL_PREFIX}
                )
find_library( GT_LOG NAMES gadgetron_toolbox_log
                NAMES_PER_DIR
                PATHS ${CMAKE_INSTALL_PREFIX}
                )

set (GT_LIBS ${GT_CPUCORE} ${GT_FFT} ${GT_NFFT} ${GT_LOG} )

if(GT_CPUCORE AND GT_FFT AND GT_NFFT AND GT_LOG)
    MESSAGE(STATUS "The Gadgetron toolboxes required for radial gridding were found. Non-cartesian encoding will be compiled.")
    set(GADGETRON_TOOLBOXES_AVAILABLE TRUE)
    set(CGADGETRON_SOURCES ${CGADGETRON_SOURCES} NonCartesianEncoding.cpp)
    add_library(cgadgetron ${CGADGETRON_SOURCES})

    if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
      add_definitions(-DGADGETRON_TOOLBOXES_AVAILABLE)
    else()
      add_compile_definitions(GADGETRON_TOOLBOXES_AVAILABLE)
    endif()
    if (Gadgetron_USE_CUDA)
        find_package(CUDA)
        if (CUDA_FOUND)
            message(STATUS "<<<<<<<<<<<<<<<<< CUDA FOUND >>>>>>>>>>>>>>>>>>>>>")
            message(STATUS "Will enable CUDA dependencies for Gadgetron.")
            target_include_directories(cgadgetron PUBLIC "${CUDA_INCLUDE_DIRS}")
	    # might need gadgetron_toolbox_gpucore as well
	    list(APPEND GT_LIBS ${CUDA_LIBRARIES}
	       ${CUDA_CUFFT_LIBRARIES}
	       ${CUDA_CUBLAS_LIBRARIES}
	       ${CUDA_CUSPARSE_LIBRARIES})
         else()
	    message(FATAL_ERROR "CUDA not found and Gadgetron is built with CUDA. Cannot compile.")
         endif()
      endif()
else()
    set(GADGETRON_TOOLBOXES_AVAILABLE FALSE)
    MESSAGE(WARNING "The Gadgetron Toolboxes required for radial gridding were NOT FOUND. Non-cartesian encoding will NOT BE compiled.")
    add_library(cgadgetron ${CGADGETRON_SOURCES})
endif()

set (cGadgetron_INCLUDE_DIR "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>$<INSTALL_INTERFACE:include>")
target_include_directories(cgadgetron PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")
# copy to parent scope
#set (cGadgetron_INCLUDE_DIR "${cGadgetron_INCLUDE_DIR}" PARENT_SCOPE)

target_include_directories(cgadgetron PUBLIC "${cGadgetron_INCLUDE_DIR}")
target_include_directories(cgadgetron PRIVATE "${FFTW3_INCLUDE_DIR}")
target_include_directories(cgadgetron PUBLIC "${ISMRMRD_INCLUDE_DIR}")

target_link_libraries(cgadgetron PUBLIC iutilities csirf)
# Add boost library dependencies
if((CMAKE_VERSION VERSION_LESS 3.5.0) OR (NOT _Boost_IMPORTED_TARGETS))
  # This is harder than it should be on older CMake versions to be able to cope with
  # spaces in filenames.
  foreach(C SYSTEM FILESYSTEM THREAD DATE_TIME CHRONO)
    target_link_libraries(cgadgetron PUBLIC optimized "${Boost_${C}_LIBRARY_RELEASE}")
    target_link_libraries(cgadgetron PUBLIC debug  "${Boost_${C}_LIBRARY_DEBUG}")
  endforeach()
else()
  # Nice and simple for recent CMake (which knows about your Boost version)
  target_link_libraries(cgadgetron PUBLIC Boost::system Boost::filesystem Boost::thread Boost::date_time Boost::chrono)
endif()

target_link_libraries(cgadgetron PUBLIC ISMRMRD::ISMRMRD)
target_link_libraries(cgadgetron PUBLIC "${FFTW3_LIBRARIES}")

if(GADGETRON_TOOLBOXES_AVAILABLE)
    target_link_libraries(cgadgetron PUBLIC ${GT_LIBS})
endif()

target_compile_definitions(cgadgetron PRIVATE _WIN32_WINNT=0x0501)

ADD_SUBDIRECTORY(tests)
