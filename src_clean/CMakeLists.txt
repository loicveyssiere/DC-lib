# ------------------------------------------------------------------------------
# Generic CMAKE file for all sub-projects
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

# Get the name of the project
get_filename_component(BENCH_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Project Name
project(${BENCH_NAME})
set(exec ${BENCH_NAME})

# Sources ----------------------------------------------------------------------
file(GLOB_RECURSE SOURCES "src/*.cc")

# Includes ---------------------------------------------------------------------
include_directories (src/headers/)
include_directories (../utils/include)
include_directories (../commons/include)
include_directories (${DC_INCLUDE})
include_directories (${HDF5_INCLUDE_DIRS})
if (${DISTRI} STREQUAL "XMPI")
    include_directories (${MPI_INCLUDE_PATH})
elseif (${DISTRI} STREQUAL "GASPI")
    include_directories (${GASPI_INCLUDE})
endif (${DISTRI} STREQUAL "XMPI")
if (${VTUNE})
    include_directories (${VTUNE_INCLUDE})
endif (${VTUNE})
if (${CILKVIEW})
    include_directories (${CILKVIEW_INCLUDE})
endif (${CILKVIEW})

# Target -----------------------------------------------------------------------
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_executable(${BENCH_NAME} ${SOURCES})

# Links ------------------------------------------------------------------------
target_link_libraries (${exec} ${DC_LIBRARIES})
target_link_libraries (${exec} ${HDF5_LIBRARIES})
target_link_libraries (${exec} utils)
target_link_libraries (${exec} commons)
if (${DISTRI} STREQUAL "XMPI")
    target_link_libraries (${exec} ${MPI_LIBRARIES})
elseif (${DISTRI} STREQUAL "GASPI")
    target_link_libraries (${exec} ${GASPI_LIBRARIES})
endif (${DISTRI} STREQUAL "XMPI")
if (${TREE})
    target_link_libraries (${exec} ${METIS_LIBRARIES})
endif (${TREE})
if (${VTUNE})
    target_link_libraries (${exec} ${VTUNE_LIBRARIES})
endif (${VTUNE})
