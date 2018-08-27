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
include_directories(${PROJECT_SOURCE_DIR}/src/Registration/cReg)
include_directories(${PROJECT_SOURCE_DIR}/src/iUtilities)

include_directories(./Testing)

set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynamicsimulation_x.cpp  Testing/tests_dynamicsimulation.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_mr_acquisition_model.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} auxiliary_input_output.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_auxiliary_testing_functions.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_auxiliary_input_output.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/auxiliary_testing_functions.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/all_simulation_tests.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} tissueparameters.cpp Testing/tests_tissueparameters.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} tissuelabelmapper.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} contrastgenerator.cpp  Testing/tests_contrastgenerator.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} phantom_input.cpp Testing/tests_phantom_input.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_encoding.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynamics.cpp Testing/tests_dynamics.cpp)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynsim_noisegenerator.cpp Testing/tests_noisegenerator.cpp)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynsim_deformer.cpp Testing/tests_dynsim_deformer.cpp)



add_library( cdynamicsimulation ${SRC_FILES_TO_COMPILE})

set (cdynamicsimulation "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")
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
  target_link_libraries(cdynamicsimulation Boost::system Boost::filesystem Boost::thread Boost::date_time Boost::chrono)
endif()

# Note: cannot use ISMRMRD_LIBRARIES on Windows as it generally contains 
# a list of filenames with spaces. There doesn't seem to be a way to pass this through (strange).
# Luckily, we know what libraries it uses
target_link_libraries(cdynamicsimulation ismrmrd)
target_link_libraries(cdynamicsimulation "${FFTW3_LIBRARIES}")
target_link_libraries(cdynamicsimulation "${HDF5_LIBRARIES}")


target_link_libraries(cdynamicsimulation cgadgetron)
target_link_libraries(cdynamicsimulation cstir)

target_link_libraries(cdynamicsimulation ${HDF5_CXX_LIBRARIES})

set( GADGETRON_LIBRARIES_TOLINK gadgetron_toolbox_cpucore ${GADGETRON_LIBRARIES} )
set( GADGETRON_LIBRARIES_TOLINK gadgetron_toolbox_cpufft ${GADGETRON_LIBRARIES_TOLINK} )
set( GADGETRON_LIBRARIES_TOLINK gadgetron_toolbox_image_analyze_io ${GADGETRON_LIBRARIES_TOLINK} )
set( GADGETRON_LIBRARIES_TOLINK gadgetron_toolbox_log ${GADGETRON_LIBRARIES_TOLINK} )

#message(STATUS "#########################GADGETRON_LIBRARIES: ${GADGETRON_LIBRARIES_TOLINK}")

target_link_libraries(cdynamicsimulation ${GADGETRON_LIBRARIES_TOLINK})

#get_property( inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#message ( STATUS "##################################### include dirs #####################################: " ${inc_dirs})

message(STATUS "XXXXXXXXXXXXXXXXXXSTIR_LIB ${STIR_LIBRARIES}")
add_executable( TestJobs main.cpp ${STIR_REGISTRIES})
target_link_libraries( TestJobs cdynamicsimulation ${STIR_LIBRARIES} )
