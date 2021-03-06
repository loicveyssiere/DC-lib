#  Copyright 2014 - UVSQ
#  Authors list: Loïc Thébault, Eric Petit
#
#  This file is part of the DC-lib.
#
#  DC-lib is free software: you can redistribute it and/or modify it under the
#  terms of the GNU Lesser General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later version.
#
#  DC-lib is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with
#  the DC-lib. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required (VERSION 2.6)

# Compilers
set (LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})
set (CMAKE_CXX_COMPILER icpc)

# Project
project (DC-lib)
if (${DISTRI} STREQUAL "XMPI")
    find_package (MPI)
endif (${DISTRI} STREQUAL "XMPI")

# Flags
set (flags "-O2 -mkl:sequential")
add_definitions (-D${DISTRI} -D${SHARED})
if (${MAX_ELEM_PER_PART} GREATER 0)
    add_definitions(-DMAX_ELEM_PER_PART=${MAX_ELEM_PER_PART})
endif (${MAX_ELEM_PER_PART} GREATER 0)
if (${BULK})
    add_definitions (-DBULK_SYNCHRONOUS)
else (${BULK})
    add_definitions (-DMULTITHREADED_COMM)
endif (${BULK})
if (${VECTO})
    add_definitions (-DDC_VEC)
endif (${VECTO})
if (${ARCHI} STREQUAL "SSE")
    set (flags "${flags} -msse")
    add_definitions (-DVEC_SIZE=2)
elseif (${ARCHI} STREQUAL "AVX")
    set (flags "${flags} -mavx")
    add_definitions (-DVEC_SIZE=4)
elseif (${ARCHI} STREQUAL "KNC")
    set (flags "${flags} -mmic")
    add_definitions (-DVEC_SIZE=8)
elseif (${ARCHI} STREQUAL "KNL")
    set (flags "${flags} -xMIC-AVX512")
    add_definitions (-DVEC_SIZE=8)
endif (${ARCHI} STREQUAL "SSE")
if (${TREE})
    add_definitions (-DTREE_CREATION)
endif (${TREE})
if (${STATS})
    add_definitions (-DSTATS)
endif (${STATS})
if (${DEBUG})
    set (flags "${flags} -g -O0")
endif (${DEBUG})
if (${VERBOSE})
    set (CMAKE_VERBOSE_MAKEFILE 1)
    set (flags "${flags} -qopt-report-phase=vec")
endif (${VERBOSE})
if (${SHARED} STREQUAL "OMP")
    set (flags "${flags} -openmp")
endif (${SHARED} STREQUAL "OMP")
set (CMAKE_CXX_FLAGS "${flags} -std=c++11")

# Headers
include_directories (../include)
include_directories (../src/headers)
include_directories (../src2/headers)
if (${DISTRI} STREQUAL "XMPI")
    include_directories (${MPI_INCLUDE_PATH})
elseif (${DISTRI} STREQUAL "GASPI")
    include_directories (${GASPI_INCLUDE})
endif (${DISTRI} STREQUAL "XMPI")
if (${TREE})
    include_directories (${METIS_INCLUDE})
endif (${TREE})

# Sources
file (GLOB_RECURSE source_list ../src*/* ../include/*)

# Library
set (lib DC)
if (${BULK})
    set (lib ${lib}_BulkSynchronous)
else (${BULK})
    set (lib ${lib}_MultithreadedComm)
endif (${BULK})
if (${VECTO})
	set (lib ${lib}_VEC_${ARCHI})
elseif (${ARCHI} STREQUAL "KNC")
	set (lib ${lib}_${ARCHI})
elseif (${ARCHI} STREQUAL "KNL")
	set (lib ${lib}_${ARCHI})
endif (${VECTO})
if (${TREE})
	set (lib ${lib}_TreeCreation)
endif (${TREE})
if (${STATS})
	set (lib ${lib}_Stats)
endif (${STATS})
if (${DEBUG})
	set (lib ${lib}_Debug)
endif (${DEBUG})
set (lib ${lib}_${DISTRI}_${SHARED})

add_library (
	${lib}
    SHARED
	${source_list}
)
