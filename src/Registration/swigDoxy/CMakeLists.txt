#========================================================================
# Author: Richard Brown
# Copyright 2016, 2017 University College London
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

# Find doxygen
FIND_PACKAGE(Doxygen)
IF(NOT DOXYGEN_FOUND)
	RETURN()
ENDIF()

# Configure and run doxygen on the C files and save log to file
CONFIGURE_FILE(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
EXECUTE_PROCESS(COMMAND ${DOXYGEN_EXECUTABLE} 
	${CMAKE_CURRENT_BINARY_DIR}/Doxyfile 
	OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.log)

# Run the doxy2swig script and save log to file
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} 
	${CMAKE_CURRENT_SOURCE_DIR}/doxy2swig.py 
	-f -t -c -a -o
	${CMAKE_CURRENT_BINARY_DIR}/xml/index.xml 
	${CMAKE_CURRENT_BINARY_DIR}/SWIG_SIRFREG_DOCUMENTATION.i
	OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/doxy2swig.log)