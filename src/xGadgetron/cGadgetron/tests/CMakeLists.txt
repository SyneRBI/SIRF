#========================================================================
# Author: Johannes Mayer
# Author: Kris Thielemans
# Author: Edoardo Pasca
# Copyright 2016 - 2020, 2022 University College London
# Copyright 2021, 2022 Science Technology Facilities Council
# Copyright 2020-2022 Physikalisch-Technische Bundesanstalt Berlin
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


function(get_zenodo_file path url md5)
        message(STATUS "Downloading (if not already done)...
                Path:     ${path}
                URL:      ${url}
                MD5:      ${md5}")
        file(DOWNLOAD "${url}" "${path}" SHOW_PROGRESS EXPECTED_MD5 "${md5}")
endfunction()

OPTION(DOWNLOAD_ZENODO_TEST_DATA "download zenodo data for tests" ON)
if (DOWNLOAD_ZENODO_TEST_DATA)
        set(MR_RPE_DATA "${CMAKE_SOURCE_DIR}/data/examples/MR/zenodo/3D_RPE_Lowres.h5")
        get_zenodo_file(
                "${MR_RPE_DATA}"
                "https://zenodo.org/record/4600937/files/3DRPE_Lowres.h5?download=1"
                "a692ba3a6c456d1ad49750717a9d5dda")
endif()

OPTION(RUN_ISMRMRD_SHEPP_LOGAN "Generate testdata compatible with ISMRMRD Version" ON)
find_program(EXE_ISMRMRD_SL NAMES ismrmrd_generate_cartesian_shepp_logan
                NAMES_PER_DIR
                PATHS ${CMAKE_INSTALL_PREFIX})

set(TESTDATA_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/data/examples/MR")
set(TESTDATA_OUTPUT_NAME_FULL "simulated_MR_2D_cartesian.h5")
set(TESTDATA_OUTPUT_NAME_GRAPPA "simulated_MR_2D_cartesian_Grappa2.h5")

add_library(MR_TESTS_CPP_AUXILIARY ${CMAKE_CURRENT_SOURCE_DIR}/mrtest_auxiliary_funs.cpp)
target_link_libraries(MR_TESTS_CPP_AUXILIARY PUBLIC csirf cgadgetron)

add_executable(MR_PROCESS_TESTDATA ${CMAKE_CURRENT_SOURCE_DIR}/mrtests_prep_testdata.cpp)
target_link_libraries(MR_PROCESS_TESTDATA PUBLIC MR_TESTS_CPP_AUXILIARY csirf cgadgetron)

if(RUN_ISMRMRD_SHEPP_LOGAN)
  if(EXE_ISMRMRD_SL)
    if (NOT EXISTS ${TESTDATA_OUTPUT_PATH})
      file(MAKE_DIRECTORY ${TESTDATA_OUTPUT_PATH})
    endif()
    set(TESTDATA_TEMP_PREFIX "testdata_tmp")
    set (TESTDATA_OUTPUT_FILES
      "${TESTDATA_OUTPUT_NAME_FULL}"
      "${TESTDATA_OUTPUT_NAME_GRAPPA}")
    set(TESTDATA_TEMP_FILES
      "${TESTDATA_TEMP_PREFIX}_1.h5"
      "${TESTDATA_TEMP_PREFIX}_2.h5")

    # We will first execute ismrmrd_generate_cartesian_shepp_logan, followed
    # by a processing step calling an executable MR_PROCESS_TESTDATA which fills a
    # default geometry orientation into the MR acquisitions to suppress warnings.

    # Sadly, this is complicated as these executables probably depend on shared
    # libraries. Their locations need to be given to the OS when running the executables.
    # This might have been done before calling CMake, but we are not sure.
    # So, we will try to determine/guess their locations.
    if (CMAKE_VERSION VERSION_LESS 3.20.0)
      message(STATUS "CMAKE_VERSION less than 3.20. Hopefully ISMRMRD shared library is in the same location as where we will install SIRF")
      # (or is already found by the OS)
      set (LIBRTPATH ${CMAKE_INSTALL_PREFIX}/lib)
    else()
      # We need ISMRMRD, Boost, HDF5 and FFTW libraries...
      # We will try to find locations based on corresponding variables.

      # separator used in path variables is different on Windows and elsewhere
      if (WIN32)
        set(PATHSEP ";")
      else()
        set(PATHSEP ":")
      endif()

      if (WIN32)
        cmake_path(REMOVE_FILENAME EXE_ISMRMRD_SL OUTPUT_VARIABLE ISMRMRDPATH)
      else()
        get_target_property(ISMRDRDCONFIGS ISMRMRD::ISMRMRD IMPORTED_CONFIGURATIONS)
        list(GET ISMRDRDCONFIGS 0 ISMRDRDCONFIG)
        get_target_property(ISMRMRDPATH ISMRMRD::ISMRMRD IMPORTED_LOCATION_${ISMRDRDCONFIG})
        cmake_path(REMOVE_FILENAME ISMRMRDPATH OUTPUT_VARIABLE ISMRMRDPATH)
        message(STATUS "ISMRMRD shared library location: ${ISMRMRDPATH}")
      endif()

      # Sadly, naming of HDF5 variables is quite unpredictable, so we try a few.
      if (HDF5_C_LIBRARY)
        cmake_path(REMOVE_FILENAME HDF5_C_LIBRARY OUTPUT_VARIABLE HDF5PATH)
      elseif (HDF5_hdf5_cpp_LIBRARY_RELEASE)
        cmake_path(REMOVE_FILENAME HDF5_hdf5_cpp_LIBRARY_RELEASE OUTPUT_VARIABLE HDF5PATH)
      endif()
      if (WIN32 AND HDF5PATH)
        # Guess the location of the HDF5 DLL from the .lib.
        set (HDF5PATH "${HDF5PATH}/../bin")
      endif()
      #message(STATUS HDF5PATH=${HDF5PATH})

      if (FFTW3_LIBRARY)
        cmake_path(REMOVE_FILENAME FFTW3_LIBRARY OUTPUT_VARIABLE FFTWPATH)
      elseif(FFTW3F_LIBRARY)
        cmake_path(REMOVE_FILENAME FFTW3F_LIBRARY OUTPUT_VARIABLE FFTWPATH)
      endif()

      set(LIBRTPATH "${CMAKE_INSTALL_PREFIX}/lib${PATHSEP}${ISMRMRDPATH}${PATHSEP}${Boost_LIBRARY_DIR_RELEASE}${PATHSEP}${FFTWPATH}${PATHSEP}${HDF5PATH}")
    endif()
    message(STATUS "path for finding shared libs/DLLs when generating MR data: ${LIBRTPATH}")

    # Now create a variable that sets the respective variable
    if (WIN32)
      # On Windows, we need to set the path.
      set (envarg "path=${LIBRTPATH};$ENV{Path}")
    elseif (APPLE)
      # On MacOS, we set DYLD_LIBRARY_PATH
      set (envarg "DYLD_LIBRARY_PATH=${LIBRTPATH}:$ENV{DYLD_LIBRARY_PATH}")
    else() 
      # On Linux, we set LD_LIBRARY_PATH
      set (envarg "LD_LIBRARY_PATH=${LIBRTPATH}:$ENV{LD_LIBRARY_PATH}")
    endif()

    # rm has been added in CMake 3.17
    if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.17")
      set(RM rm)
    else()
      set(RM remove)
    endif()

    # Finally, add the custom target that will create the files
    add_custom_target(GENERATE_SIMULATED_TESTDATA ALL
      COMMAND ${CMAKE_COMMAND} -E ${RM} -f ${TESTDATA_OUTPUT_FILES} ${TESTDATA_TEMP_FILES}
      COMMAND ${CMAKE_COMMAND} -E env "${envarg}" ${EXE_ISMRMRD_SL} -r 2 -o "${TESTDATA_TEMP_PREFIX}_1.h5"
      COMMAND ${CMAKE_COMMAND} -E env "${envarg}" ${EXE_ISMRMRD_SL} -w 64 -a 2 -o "${TESTDATA_TEMP_PREFIX}_2.h5"
      COMMAND ${CMAKE_COMMAND} -E env "${envarg}" $<TARGET_FILE:MR_PROCESS_TESTDATA> "${TESTDATA_TEMP_PREFIX}_1.h5" ${TESTDATA_OUTPUT_NAME_FULL}
      COMMAND ${CMAKE_COMMAND} -E env "${envarg}" $<TARGET_FILE:MR_PROCESS_TESTDATA> "${TESTDATA_TEMP_PREFIX}_2.h5" ${TESTDATA_OUTPUT_NAME_GRAPPA}
      COMMAND ${CMAKE_COMMAND} -E ${RM} -f ${TESTDATA_TEMP_FILES}
      WORKING_DIRECTORY ${TESTDATA_OUTPUT_PATH})

    # adding a dependency of PROCESS_SIMULATED_TESTDATA on MR_PROCESS_TESTDATA ensures the
    # executable is built before it is called by the custom target.
    add_dependencies(GENERATE_SIMULATED_TESTDATA MR_PROCESS_TESTDATA)
   
  else()
    message(FATAL_ERROR "We could not find the ismrmrd_generate_cartesian_shepp_logan executable. Either set the corresponding variable EXE_ISMRMRD_SL=<path_to>ismrmrd_generate_cartesian_shepp_logan_executable, or disable the generation of the Shepp Logan phantom by setting RUN_ISMRMRD_SHEPP_LOGAN=OFF" )
  endif()
else()
    # do nothing but add target as it will be a dependency later
    add_custom_target(GENERATE_SIMULATED_TESTDATA ALL COMMAND "")
  message(WARNING "Will not run the Shepp Logan Generation." )
endif()


add_executable(MR_TESTS_CPLUSPLUS ${CMAKE_CURRENT_SOURCE_DIR}/mrtests.cpp)
add_dependencies(MR_TESTS_CPLUSPLUS GENERATE_SIMULATED_TESTDATA)

target_link_libraries(MR_TESTS_CPLUSPLUS PUBLIC MR_TESTS_CPP_AUXILIARY csirf cgadgetron)
INSTALL(TARGETS MR_TESTS_CPLUSPLUS DESTINATION bin)


# generate filenames that are passed to the C++ test exe and install testdata into the shared directory.
set(MR_TESTDATA_PATH "${SHARE_DIR}/data/examples/MR")
if(RUN_ISMRMRD_SHEPP_LOGAN)
  install(FILES  "${TESTDATA_OUTPUT_PATH}/${TESTDATA_OUTPUT_NAME_FULL}"
                 "${TESTDATA_OUTPUT_PATH}/${TESTDATA_OUTPUT_NAME_GRAPPA}"
          DESTINATION "${MR_TESTDATA_PATH}")
  set(FILENAME_CARTESIAN_TESTDATA "${MR_TESTDATA_PATH}/${TESTDATA_OUTPUT_NAME_FULL}")
else()
    set(FILENAME_CARTESIAN_TESTDATA "${MR_TESTDATA_PATH}/simulated_MR_2D_cartesian.h5")
endif()
if(GADGETRON_TOOLBOXES_AVAILABLE)
  set(FILENAME_RPE_TESTDATA "${MR_TESTDATA_PATH}/zenodo/3D_RPE_Lowres.h5")
else()
  set(FILENAME_RPE_TESTDATA "")
endif()

ADD_TEST(NAME MR_TESTS_CPLUSPLUS
         COMMAND MR_TESTS_CPLUSPLUS "${FILENAME_CARTESIAN_TESTDATA}" "${FILENAME_RPE_TESTDATA}"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

