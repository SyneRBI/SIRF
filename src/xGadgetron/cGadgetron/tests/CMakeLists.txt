#========================================================================
# Author: Johannes Mayer
# Copyright 2016 - 2020 University College London
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

set(TESTDATA_TEMP_PREFIX "testdata_tmp")

if(RUN_ISMRMRD_SHEPP_LOGAN)
  if(EXE_ISMRMRD_SL)
    if (NOT EXISTS ${TESTDATA_OUTPUT_PATH})
      file(MAKE_DIRECTORY ${TESTDATA_OUTPUT_PATH})
    endif()
    # rm has been added in CMake 3.17 
    if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.17")
      add_custom_target(GENERATE_SIMULATED_TESTDATA ALL
          COMMAND ${CMAKE_COMMAND} -E rm -f "${TESTDATA_OUTPUT_NAME_FULL}"
          COMMAND ${CMAKE_COMMAND} -E rm -f "${TESTDATA_OUTPUT_NAME_GRAPPA}"
          COMMAND ${CMAKE_COMMAND} -E rm -f "${TESTDATA_TEMP_PREFIX}_1.h5"
          COMMAND ${CMAKE_COMMAND} -E rm -f "${TESTDATA_TEMP_PREFIX}_2.h5"
          COMMAND ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH} ${EXE_ISMRMRD_SL} -r 2 -o "${TESTDATA_TEMP_PREFIX}_1.h5"
          COMMAND ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH} ${EXE_ISMRMRD_SL} -w 64 -a 2 -o "${TESTDATA_TEMP_PREFIX}_2.h5"
          WORKING_DIRECTORY ${TESTDATA_OUTPUT_PATH})
    else()
      add_custom_target(GENERATE_SIMULATED_TESTDATA ALL
          COMMAND ${CMAKE_COMMAND} -E remove -f "${TESTDATA_OUTPUT_NAME_FULL}"
          COMMAND ${CMAKE_COMMAND} -E remove -f "${TESTDATA_OUTPUT_NAME_GRAPPA}"
          COMMAND ${CMAKE_COMMAND} -E remove -f "${TESTDATA_TEMP_PREFIX}_1.h5"
          COMMAND ${CMAKE_COMMAND} -E remove -f "${TESTDATA_TEMP_PREFIX}_2.h5"
          COMMAND ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH} ${EXE_ISMRMRD_SL} -r 2 -o "${TESTDATA_TEMP_PREFIX}_1.h5"
          COMMAND ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH} ${EXE_ISMRMRD_SL} -w 64 -a 2 -o "${TESTDATA_TEMP_PREFIX}_2.h5"
          WORKING_DIRECTORY ${TESTDATA_OUTPUT_PATH})
    endif()
    
  else()
    message(FATAL_ERROR "We could not find the ismrmrd_generate_cartesian_shepp_logan executable. Either set the corresponding variable EXE_ISMRMRD_SL=<path_to>ismrmrd_generate_cartesian_shepp_logan_executable, or disable the generation of the Shepp Logan phantom by setting RUN_ISMRMRD_SHEPP_LOGAN=OFF" )
  endif()
else()
    # do nothing but add target as it will be a dependence later
    add_custom_target(GENERATE_SIMULATED_TESTDATA ALL COMMAND "")
  message(WARNING "Will not run the Shepp Logan Generation." )
endif()

add_library(MR_TESTS_CPP_AUXILIARY ${CMAKE_CURRENT_SOURCE_DIR}/mrtest_auxiliary_funs.cpp)
target_link_libraries(MR_TESTS_CPP_AUXILIARY PUBLIC csirf cgadgetron)

add_executable(MR_PROCESS_TESTDATA ${CMAKE_CURRENT_SOURCE_DIR}/mrtests_prep_testdata.cpp)
add_dependencies(MR_PROCESS_TESTDATA GENERATE_SIMULATED_TESTDATA)
target_link_libraries(MR_PROCESS_TESTDATA PUBLIC MR_TESTS_CPP_AUXILIARY csirf cgadgetron)
set(MR_PROCESS_TESTDATA_EXE "$<TARGET_FILE:MR_PROCESS_TESTDATA>")


if(RUN_ISMRMRD_SHEPP_LOGAN)
  # this pre-processing step calls an executable MR_PROCESS_TESTDATA which fills a
  # default geometry orientation into the MR acquisitions to suppress warnings.
  add_custom_target(PROCESS_SIMULATED_TESTDATA ALL
          COMMAND ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH} ${MR_PROCESS_TESTDATA_EXE} "${TESTDATA_TEMP_PREFIX}_1.h5" ${TESTDATA_OUTPUT_NAME_FULL}
          COMMAND ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH} ${MR_PROCESS_TESTDATA_EXE} "${TESTDATA_TEMP_PREFIX}_2.h5" ${TESTDATA_OUTPUT_NAME_GRAPPA}
          COMMAND ${CMAKE_COMMAND} -E rm -f "${TESTDATA_TEMP_PREFIX}_1.h5"
          COMMAND ${CMAKE_COMMAND} -E rm -f "${TESTDATA_TEMP_PREFIX}_2.h5"
          WORKING_DIRECTORY ${TESTDATA_OUTPUT_PATH})
else()
        add_custom_target(PROCESS_SIMULATED_TESTDATA ALL
        COMMAND "")
endif()

# adding a dependency of PROCESS_SIMULATED_TESTDATA on MR_PROCESS_TESTDATA ensures the
# executable is built before it is called by the custom target.
add_dependencies(PROCESS_SIMULATED_TESTDATA MR_PROCESS_TESTDATA)

add_executable(MR_TESTS_CPLUSPLUS ${CMAKE_CURRENT_SOURCE_DIR}/mrtests.cpp)
add_dependencies(MR_TESTS_CPLUSPLUS GENERATE_SIMULATED_TESTDATA PROCESS_SIMULATED_TESTDATA)

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

