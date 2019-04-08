#========================================================================
# Author: Richard Brown
# Copyright 2018 - 2019 University College London
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


FUNCTION(read_methods file, methods_list, NR_VERSION, NR_VERSION_to_use)
  # Read and break into a list separated by new lines
  FILE(READ ${methods_file} contents)
  # Make file into list (new elements at newline)
  STRING(REGEX REPLACE ";" "\\\\;" contents "${contents}")
  STRING(REGEX REPLACE "\n" ";" contents "${contents}")
  # loop over all lines, ignore those starting with '#'
  FOREACH(elem ${contents})
    IF(${elem} MATCHES "#")
      CONTINUE()
    ENDIF()

    # Check against known versions
    IF(${elem} MATCHES "Known versions")
      IF(${elem} MATCHES ${NR_VERSION})
        SET("NR_VERSION_to_use" "${NR_VERSION}" PARENT_SCOPE)
      ELSE()
        STRING(REGEX MATCHALL "[^,]+" "version_list" "${elem}")
        LIST(GET "version_list" "-1" NEW_NR_VERSION)
        MESSAGE("version of NR not found (${NR_VERSION}), using most recent methods (${NEW_NR_VERSION})")
        SET("NR_VERSION_to_use" "${NEW_NR_VERSION}" PARENT_SCOPE)
      ENDIF()
    # The rest are the table of methods, add to list
    ELSE()
      LIST(APPEND "table" ${elem})
    ENDIF()
  ENDFOREACH(elem ${REG_executables})
  SET("methods_list" "${table}" PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(get_parser_methods methods_file class_name NR_VERSION result_variable)
  read_methods(methods_file, methods_list, ${NR_VERSION}, NR_VERSION_to_use)
  FOREACH(method_line ${methods_list})  
    # make into list
    IF("${method_line}" MATCHES "${NR_VERSION_to_use}")
      STRING(REGEX MATCHALL "[^,]+" "elem" "${method_line}")
      LIST(GET "elem" 0 "command")
      STRING(APPEND output "    parser.add_key(\"${command}\",&${class_name}<dataType>::${command});\n")
    ENDIF()
  ENDFOREACH()
  SET("${result_variable}" "${output}" PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(get_method_type_as_string type arg_num str)
  STRING(REPLACE " " "" type "${type}")
  STRING(REPLACE "," "" type "${type}")
  IF("${type}" MATCHES "-")
    return()
  ELSEIF("${type}" MATCHES "str")
    SET(stri "arg${arg_num}.c_str()")
  ELSEIF("${type}" MATCHES "int")
    SET(stri "stoi(arg${arg_num})")
  ELSEIF("${type}" MATCHES "float")
    SET(stri "stof(arg${arg_num})")
  ELSEIF("${type}" MATCHES "double")
    SET(stri "stod(arg${arg_num})")
  ELSEIF("${type}" MATCHES "unsigned")
    SET(stri "unsigned(stoi(arg${arg_num}))")
  ELSEIF("${type}" MATCHES "dataType")
    SET(stri "dataType(stod(arg${arg_num}))")
  ELSEIF("${type}" MATCHES "bool")
    SET(stri "bool(stoi(arg${arg_num}))")
  ENDIF()
  IF(${arg_num} MATCHES "1")
    SET("${str}" "${stri}" PARENT_SCOPE)
  ELSE()
    SET("${str}" ", ${stri}" PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

FUNCTION(get_runtime_methods methods_file NR_VERSION result_variable)
  read_methods(methods_file, methods_list, ${NR_VERSION}, NR_VERSION_to_use)
  FOREACH(method_line ${methods_list})  
    # make into list
    IF("${method_line}" MATCHES "${NR_VERSION_to_use}")
      STRING(REGEX MATCHALL "[^,]+" "elem" "${method_line}")
      LIST(GET "elem" 0 "command")
      LIST(GET "elem" 1 "type1")
      LIST(GET "elem" 2 "type2")
      IF(NOT output)
        STRING(APPEND output "        if     ")
      ELSE()
        STRING(APPEND output "        else if")
      ENDIF()
      STRING(APPEND output " (strcmp(par.c_str(),\"${command}\")== 0) _registration_sptr->${command}(")
      get_method_type_as_string("${type1}" "1" "type1_str")
      get_method_type_as_string("${type2}" "2" "type2_str")
      STRING(APPEND output "${type1_str}${type2_str}")
      STRING(APPEND output ");\n")

    ENDIF()
  ENDFOREACH()
  SET("${result_variable}" "${output}" PARENT_SCOPE)
ENDFUNCTION()