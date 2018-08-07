#========================================================================
# Author: Johannes Mayer
# Date: 15.03.2018
#=========================================================================

set(CMAKE_POSITION_INDEPENDENT_CODE True)

if (CMAKE_VERSION VERSION_LESS "3.1")
  set (CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
  else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

find_package(ISMRMRD)
if (ISMRMRD_FOUND)
  set(CMAKE_MODULE_PATH "${ISMRMRD_DIR}")

  find_package(FFTW3 COMPONENTS single REQUIRED)

  if (SIRF_INSTALL_DEPENDENCIES AND WIN32)
    set(FFTW3_DLL_DIR ${FFTW3_INCLUDE_DIR})
    message(STATUS "Install FFTW3 libraries from ${FFTW3_DLL_DIR} ")
    file(GLOB FFTW_DLLS "${FFTW3_DLL_DIR}/*.dll")
    install( FILES ${FFTW_DLLS} DESTINATION bin )
  endif()

  include_directories("${ISMRMRD_INCLUDE_DIR}")
  link_directories("${ISMRMRD_LIBRARY_DIRS}")

endif()

######### Find The Gadgetron ############
find_path(GADGETRON_CMAKE_MODULES FindGadgetron.cmake HINTS $ENV{GADGETRON_HOME}/share/gadgetron/cmake) 
if (NOT GADGETRON_CMAKE_MODULES)
  MESSAGE(FATAL_ERROR "GADGETRON_CMAKE_MODULES cannot be found. 
   Try to set GADGETRON_HOME environment variable.")
endif(NOT GADGETRON_CMAKE_MODULES)

set(CMAKE_MODULE_PATH "${GADGETRON_CMAKE_MODULES}")

## find the package
find_package(Gadgetron REQUIRED)
include_directories(${GADGETRON_INCLUDE_DIR})
link_directories(${GADGETRON_LIB_DIR})


#########################################

find_package(HDF5 1.8 COMPONENTS C CXX HL REQUIRED)


include_directories(${PROJECT_SOURCE_DIR}/src/common/include)
find_package(STIR REQUIRED)
ADD_SUBDIRECTORY(cDynamicSimulation)

set(cDynamicSimulation_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cDynamicSimulation")

