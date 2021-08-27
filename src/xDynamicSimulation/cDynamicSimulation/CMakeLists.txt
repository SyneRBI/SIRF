#========================================================================
# Author: Johannes Mayer
# Date: 15.03.2018
#=========================================================================

add_definitions(-DBOOST_ALL_NO_LIB)

if (SIRF_INSTALL_DEPENDENCIES AND WIN32)
    set(Boost_DLL_DIR ${Boost_LIBRARY_DIR_RELEASE})
    message(STATUS "Install boost shared libraries from ${Boost_DLL_DIR} ")
	foreach (__boost_lib system filesystem thread date_time chrono)
	   file(GLOB Boost_DLL "${Boost_DLL_DIR}/boost_${__boost_lib}*.dll")
       install( FILES ${Boost_DLL} DESTINATION bin )
	endforeach()
  endif()

include_directories(${PROJECT_SOURCE_DIR}/src/common/include)
include_directories(${PROJECT_SOURCE_DIR}/src/xGadgetron/cGadgetron)
include_directories(${PROJECT_SOURCE_DIR}/src/xSTIR/cSTIR)
include_directories(${PROJECT_SOURCE_DIR}/src/Registration/Reg)
include_directories(${PROJECT_SOURCE_DIR}/src/iUtilities)
include_directories("${GADGETRON_INCLUDE_DIR}/../")


set(SRC_FILES_TO_COMPILE 
  cdynamicsimulation.cpp
  dynamicsimulation_x.cpp  
  auxiliary_input_output.cpp 
  tissueparameters.cpp
  tissuelabelmapper.cpp 
  contrastgenerator.cpp 
  phantom_input.cpp 
  dynamics.cpp 
  dynsim_noisegenerator.cpp 
  dynsim_deformer.cpp 
  volume_orientator.cpp)

add_library( cdynamicsimulation ${SRC_FILES_TO_COMPILE})
set (cdynamicsimulation "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")

TARGET_INCLUDE_DIRECTORIES(cdynamicsimulation PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>$<INSTALL_INTERFACE:include/>"
)

target_include_directories(cdynamicsimulation PUBLIC "${HDF5_INCLUDE_DIRS}")

target_link_libraries(cdynamicsimulation ISMRMRD::ISMRMRD iutilities cgadgetron cstir Reg "${HDF5_LIBRARIES}" ${HDF5_CXX_LIBRARIES})

# Add boost library dependencies
if((CMAKE_VERSION VERSION_LESS 3.5.0) OR (NOT _Boost_IMPORTED_TARGETS))
  # This is harder than it should be on older CMake versions to be able to cope with
  # spaces in filenames.
  foreach(C SYSTEM FILESYSTEM THREAD DATE_TIME CHRONO)
    target_link_libraries(cdynamicsimulation optimized "${Boost_${C}_LIBRARY_RELEASE}")
    target_link_libraries(cdynamicsimulation debug  "${Boost_${C}_LIBRARY_DEBUG}")
  endforeach()
else()
  # Nice and simple for recent CMake (which knows about your Boost version)
  target_link_libraries(cdynamicsimulation Boost::system Boost::filesystem Boost::thread Boost::date_time Boost::chrono)
endif()

# Note: cannot use ISMRMRD_LIBRARIES on Windows as it generally contains 
# a list of filenames with spaces. There doesn't seem to be a way to pass this through (strange).
# Luckily, we know what libraries it uses

add_subdirectory(tests)
