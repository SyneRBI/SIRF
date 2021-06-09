#========================================================================
# Author: Johannes Mayer
# Copyright 2016 - 2020 University College London
# Copyright 2020 Physikalisch-Technische Bundesanstalt Berlin
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

OPTION(RERUN_ISMRMRD_SHEPP_LOGAN "Generate testdata compatible with ISMRMRD Version" ON)
find_program(EXE_ISMRMRD_SL NAMES ismrmrd_generate_cartesian_shepp_logan
                NAMES_PER_DIR
                PATHS ${CMAKE_INSTALL_PREFIX})


set(TESTDATA_OUTPUT_PATH "${SIRF_SOURCE_DIR}/data/examples/MR/")
set(TESTDATA_OUTPUT_NAME_FULL "simulated_MR_2D_cartesian.h5")
set(TESTDATA_OUTPUT_NAME_GRAPPA "simulated_MR_2D_cartesian_Grappa2.h5")

set(TESTDATA_TEMP_PREFIX "testdata_tmp")

if(EXE_ISMRMRD_SL AND RERUN_ISMRMRD_SHEPP_LOGAN)
    add_custom_target(GENERATE_SIMULATED_TESTDATA ALL
        COMMAND rm -f ${TESTDATA_OUTPUT_NAME_FULL}
        COMMAND rm -f ${TESTDATA_OUTPUT_NAME_GRAPPA}
        COMMAND env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$LD_LIBRARY_PATH ${EXE_ISMRMRD_SL} -r 2 -o "${TESTDATA_TEMP_PREFIX}_1.h5"
        COMMAND env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$LD_LIBRARY_PATH ${EXE_ISMRMRD_SL} -w 64 -a 2 -o "${TESTDATA_TEMP_PREFIX}_2.h5"
        WORKING_DIRECTORY ${TESTDATA_OUTPUT_PATH})
elseif(RERUN_ISMRMRD_SHEPP_LOGAN)
    message(ERROR "We could not find the ISMRMRD executable for the Shepp Logan Generation." )
else()
    add_custom_target(GENERATE_SIMULATED_TESTDATA ALL
        COMMAND "")
endif()

add_library(MR_TESTS_CPP_AUXILIARY ${CMAKE_CURRENT_SOURCE_DIR}/mrtest_auxiliary_funs.cpp)
target_link_libraries(MR_TESTS_CPP_AUXILIARY PUBLIC csirf cgadgetron)

add_executable(MR_PROCESS_TESTDATA ${CMAKE_CURRENT_SOURCE_DIR}/mrtests_prep_testdata.cpp)
add_dependencies(MR_PROCESS_TESTDATA GENERATE_SIMULATED_TESTDATA)
target_link_libraries(MR_PROCESS_TESTDATA PUBLIC MR_TESTS_CPP_AUXILIARY csirf cgadgetron)

if(RERUN_ISMRMRD_SHEPP_LOGAN)
        add_custom_target(PROCESS_SIMULATED_TESTDATA ALL
                COMMAND MR_PROCESS_TESTDATA "${TESTDATA_TEMP_PREFIX}_1.h5" ${TESTDATA_OUTPUT_NAME_FULL}
                COMMAND MR_PROCESS_TESTDATA "${TESTDATA_TEMP_PREFIX}_2.h5" ${TESTDATA_OUTPUT_NAME_GRAPPA}
                COMMAND rm -f "${TESTDATA_TEMP_PREFIX}_1.h5"
                COMMAND rm -f "${TESTDATA_TEMP_PREFIX}_2.h5"
                WORKING_DIRECTORY ${TESTDATA_OUTPUT_PATH})
else()
        add_custom_target(PROCESS_SIMULATED_TESTDATA ALL
        COMMAND "")
endif()
add_dependencies(PROCESS_SIMULATED_TESTDATA MR_PROCESS_TESTDATA)

add_executable(MR_TESTS_CPLUSPLUS ${CMAKE_CURRENT_SOURCE_DIR}/mrtests.cpp)
add_dependencies(MR_TESTS_CPLUSPLUS GENERATE_SIMULATED_TESTDATA PROCESS_SIMULATED_TESTDATA)

target_link_libraries(MR_TESTS_CPLUSPLUS PUBLIC MR_TESTS_CPP_AUXILIARY csirf cgadgetron)
INSTALL(TARGETS MR_TESTS_CPLUSPLUS DESTINATION bin)

ADD_TEST(NAME MR_TESTS_CPLUSPLUS
         COMMAND MR_TESTS_CPLUSPLUS
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

