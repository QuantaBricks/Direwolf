
# Minimal CMake build for mctc-gcp (LGPL-3.0-or-later, grimme-lab/gcp).
# The upstream project ships a meson.build, but this tree has no meson/ninja
# installed and every other third_party/ dependency here (simple-dftd3,
# mctc-lib, toml-f, libxc) already builds via CMake - this file just
# compiles the same 4 source files meson's src/meson.build lists, as a
# shared library, matching simple-dftd3's own CMakeLists.txt pattern
# (BUILD_SHARED_LIBS=ON, dynamic link so Engine's own code stays out of
# gcp's LGPL boundary - see third_party/simple-dftd3/build for the sibling
# example this mirrors).
cmake_minimum_required(VERSION 3.14)
project(gcp LANGUAGES Fortran VERSION 2.3.2)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_library(gcp
  src/gcp_precision.f90
  src/gcp_strings.f90
  src/gcp_version.f90
  src/gcp.f90
)
set_target_properties(gcp PROPERTIES
  OUTPUT_NAME "gcp"
  VERSION ${PROJECT_VERSION}
  SOVERSION 2
)
target_include_directories(gcp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
)

install(TARGETS gcp DESTINATION lib)
