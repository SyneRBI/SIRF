#========================================================================
# Author: Kris Thielemans
# Author: Jannis Fisher
# Copyright 2016 ETH Zurich
# Copyright 2017, 2018, 2021 University College London
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


# This file is largely based on the equivalent code in  STIR.

## Build doxygen documentation
# Has to be after previous things such that all variables are filled
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})
if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()
endif()

if(DOXYGEN_FOUND)
    if(DOXYGEN_DOT_FOUND)
        set(HAVE_DOT 1)
        message(STATUS "Graphviz dot found. It will be used for the Doxygen documentation.")
    else()
        set(HAVE_DOT 0)
        message(WARNING "Graphviz dot is not found. It is recommended to install Graphviz for nicer Doxygen documentation")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # attempt to find include directories, but this doesn't work for SIRF
    # as they are in xSTIR etc etc.
    #get_property(DOX_STRIP_FROM_INC_PATH DIRECTORY ${CMAKE_PROJECT_SOURCE_DIR}/src PROPERTY INCLUDE_DIRECTORIES)
    #string(REPLACE ";" " \\\n\t\t\t\t" DOX_STRIP_FROM_INC_PATH "${DOX_STRIP_FROM_INC_PATH}")

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    # The following lines are modifications of Jannis Fischer's original. The modifications were inspired by
    # https://samthursfield.wordpress.com/2015/11/21/cmake-dependencies-between-targets-and-files-and-custom-commands/

    # Add a custom command that does everything, but also generates a doxygen.stamp file.
    # This should prevent doxygen being re-run when not needed.
    add_custom_command(
        OUTPUT doxygen.stamp
        DEPENDS ${doxyfile} #${ALL_HEADERS} TODO
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        COMMAND ${CMAKE_COMMAND} -E touch doxygen.stamp
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
    # Now add a custom target that depends on the doxygen.stamp.
    # If BUILD_DOCUMENTATION, we add this target to the dependency list of "ALL"
    # to make sure it always gets run.
    # Otheriwse, the user has to build the target manually.
    if (BUILD_DOCUMENTATION)
      set(DOXY_ALL ALL)
    endif()

    add_custom_target(
        RUN_DOXYGEN ${DOXY_ALL}
        DEPENDS doxygen.stamp)

    if (BUILD_DOCUMENTATION)
      # Doxyfile.in currently tells doxygen to output in PROJECT_BINARY_DIR
      install(DIRECTORY ${PROJECT_BINARY_DIR}/html/
        COMPONENT DOC
        DESTINATION ${SHARE_DIR}/doc/doxygen)
    else()
      # CMake currently doesn't allow adding a dependency on a custom target.
      # As we don't know if the user has executaged the `RUN_DOXYGEN` target,
      # we hide this install a bit by adding it do a new component.
      # There has to be a better way...
      install(DIRECTORY ${PROJECT_BINARY_DIR}/html/
        COMPONENT DOXYGEN EXCLUDE_FROM_ALL
        DESTINATION ${SHARE_DIR}/doc/doxygen)
    endif()
        
endif() # DOXYGEN_FOUND
