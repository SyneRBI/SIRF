#========================================================================
# Author: Johannes Mayer
# Date: 15.03.2018
#=========================================================================


include_directories(${PROJECT_SOURCE_DIR}/src/xGadgetron/cGadgetron)

get_property( inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message ( STATUS "##################################### include dirs #####################################: " ${inc_dirs})

add_library(cdynamicSimulation  dynamicsimulation_x.cpp dynamicsimulation_x.h)