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

#########################################

# find_package(HDF5 1.8 COMPONENTS C CXX HL REQUIRED)


# include_directories(${PROJECT_SOURCE_DIR}/src/common/include)
# find_package(STIR REQUIRED)
add_subdirectory(cDynamicSimulation)

#if (PYTHON_FOUND AND BUILD_PYTHON)
#  add_subdirectory(pDynamicSimulationon)
#endif()

# set(cDynamicSimulation_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cDynamicSimulation")

