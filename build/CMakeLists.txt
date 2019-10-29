#  Copyright 2019 - UVSQ
#  Authors list: Nathalie Möller, Eric Petit
#
#  This file is part of the FMM-viz.
#
#  FMM-viz is free software: you can redistribute it and/or modify it under the
#  terms of the GNU Lesser General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later version.
#
#  FMM-viz is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with
#  the FMM-viz. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required (VERSION 2.6)

PROJECT (FMM-viz)

# Gaspi Lib
set (GASPI_VERSION "gpi")
set (GASPI_DIR "/home/nath/opt/${GASPI_VERSION}")
set (GASPI_LIBRARIES "-lGPI2-dbg -lpthread -L${GASPI_DIR}/lib64")

# FMM lib
set (FMM_LB_INCLUDE_DIR "../../FMM-lib/src/LB")
set (FMM_lib /path/to/lib/FMM-lib/lib/libfmm.so)

# OpenGL
set (OPENGL_LIBRAIRIES "-lGL -lGLU -lglut -lm")

# config
set (BIN_DIR "../bin")
set (CMAKE_CXX_COMPILER g++)

# sources
file (GLOB_RECURSE VIZ_src_list ../src/*)

# set compile FLAGS
set (CMAKE_CXX_FLAGS "-no-pie -std=c++11 -Wall -fopenmp")
set (MPI_COMPILE_FLAGS "-no-pie -std=c++11 -Wall -fopenmp")

# set executables
add_executable (viz  ${VIZ_src_list})

target_include_directories(viz PUBLIC 
	"../src/"
	${FMM_LB_INCLUDE_DIR}
	"${GASPI_DIR}/include"
)

target_link_libraries(viz ${FMM_lib})

# set bin repository
set_target_properties(viz
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)

# verbosity
set (CMAKE_VERBOSE_MAKEFILE 0)

# MPI
find_package (MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

target_link_libraries(viz 
	${MPI_LIBRARIES} 
	${OPENGL_LIBRAIRIES}
	${GASPI_LIBRARIES}
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(viz
  PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(viz
  PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
