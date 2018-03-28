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
include_directories(${PROJECT_SOURCE_DIR}/src/iUtilities)


set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynamicsimulation_x.cpp dynamicsimulation_x.h)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/auxiliary_testing_functions.cpp Testing/auxiliary_testing_functions.h)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/dynamicsimulation_tests.cpp Testing/dynamicsimulation_tests.h)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} tissueparameters.cpp tissueparameters.h Testing/tests_tissueparameters.cpp Testing/tests_tissueparameters.h)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} tissuelabelmapper.cpp tissuelabelmapper.h)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} contrastgenerator.cpp contrastgenerator.h Testing/tests_contrastgenerator.cpp Testing/tests_contrastgenerator.h)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} phantom_input.cpp phantom_input.h Testing/tests_phantom_input.cpp Testing/tests_phantom_input.h)



add_library( cdynamicsimulation ${SRC_FILES_TO_COMPILE})

set (dynamicsimulation "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")
# copy to parent scope

target_include_directories(cdynamicsimulation PUBLIC "${cDynamicSimulation_INCLUDE_DIR}")
target_include_directories(cdynamicsimulation PRIVATE "${FFTW3_INCLUDE_DIR}")
target_include_directories(cdynamicsimulation PUBLIC "${HDF5_INCLUDE_DIRS}")

target_link_libraries(cdynamicsimulation iutilities)
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
  target_link_libraries(dynamicsimulation Boost::system Boost::filesystem Boost::thread Boost::date_time Boost::chrono)
endif()

# Note: cannot use ISMRMRD_LIBRARIES on Windows as it generally contains 
# a list of filenames with spaces. There doesn't seem to be a way to pass this through (strange).
# Luckily, we know what libraries it uses
target_link_libraries(cdynamicsimulation ismrmrd)
target_link_libraries(cdynamicsimulation "${FFTW3_LIBRARIES}")
target_link_libraries(cdynamicsimulation "${HDF5_LIBRARIES}")
target_link_libraries(cdynamicsimulation ${HDF5_CXX_LIBRARIES})


#get_property( inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#message ( STATUS "##################################### include dirs #####################################: " ${inc_dirs})


add_executable( TestJobs main.cpp)
target_link_libraries( TestJobs cdynamicsimulation)
