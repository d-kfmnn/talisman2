cmake_minimum_required(VERSION 3.16)
project(talisman C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

# ============================================================
# Dependencies via FetchContent
# ============================================================

# FLINT
include(ExternalProject)

ExternalProject_Add(flint_external
  GIT_REPOSITORY https://github.com/wbhart/flint2.git
  GIT_TAG v3.2.0
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ./bootstrap.sh && ./configure
  BUILD_COMMAND make -j
  INSTALL_COMMAND ""
  BUILD_BYPRODUCTS
    ${CMAKE_BINARY_DIR}/flint_external-prefix/src/flint_external/libflint.so
)


# Correct path to produced shared library
set(FLINT_LIB
  ${CMAKE_BINARY_DIR}/flint_external-prefix/src/flint_external/libflint.so
)

add_library(flint SHARED IMPORTED GLOBAL)

set_target_properties(flint PROPERTIES
  IMPORTED_LOCATION ${FLINT_LIB}
)

add_dependencies(flint flint_external)



include(ExternalProject)

ExternalProject_Add(abc_external
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/includes/abc
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/includes/abc/rebuild_without_cadical.sh
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_BYPRODUCTS
    ${CMAKE_CURRENT_SOURCE_DIR}/includes/abc/libabc.a
)


add_library(abc STATIC IMPORTED GLOBAL)

set_target_properties(abc PROPERTIES
  IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/includes/abc/libabc.a
)

add_dependencies(abc abc_external)





# CaDiCaL
set(CADICAL_SRC_DIR ${CMAKE_BINARY_DIR}/cadical-src)

ExternalProject_Add(cadical_external
  GIT_REPOSITORY https://github.com/arminbiere/cadical.git
  GIT_TAG master
  BUILD_IN_SOURCE 1

  # ---- PATCH STEP ----
  PATCH_COMMAND
    ${CMAKE_COMMAND} -E echo "Patching CaDiCaL options.hpp" &&
    ${CMAKE_COMMAND} -E sleep 1 &&
    sed -i
    "s/OPTION( factor,            1,  0,  1,0,1,1, \"bounded variable addition\")/OPTION( factor,            0,  0,  1,0,1,1, \"bounded variable addition\")/"
    src/options.hpp

  CONFIGURE_COMMAND ./configure
  BUILD_COMMAND make -j
  INSTALL_COMMAND ""
  BUILD_BYPRODUCTS
    ${CMAKE_BINARY_DIR}/cadical_external-prefix/src/cadical_external/build/libcadical.a
)


# Create imported target manually
add_library(cadical STATIC IMPORTED GLOBAL)

set_target_properties(cadical PROPERTIES
  IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/cadical_external-prefix/src/cadical_external/build/libcadical.a
)

add_dependencies(cadical cadical_external)





# PBLib
FetchContent_Declare(
  pblib
  GIT_REPOSITORY https://github.com/master-keying/pblib.git
  GIT_TAG master
)
FetchContent_MakeAvailable(pblib)

# cmsgen
FetchContent_Declare(
  cmsgen
  GIT_REPOSITORY https://github.com/meelgroup/cmsgen.git
  GIT_TAG master
)
FetchContent_MakeAvailable(cmsgen)

# ============================================================
# Your Sources
# ============================================================

file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS src/*.cpp)
list(REMOVE_ITEM SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/cadical_tester.cpp")

add_library(aiger_obj OBJECT includes/aiger/aiger.c)
set_source_files_properties(includes/aiger/aiger.c PROPERTIES LANGUAGE C)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(aiger_obj PRIVATE -w)
endif()

add_executable(talisman
    ${SRC_FILES}
    $<TARGET_OBJECTS:aiger_obj>
)

# Include directory added manually to your executable
target_include_directories(talisman PRIVATE
  ${FLINT_INSTALL_DIR}/include
  ${cmsgen_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/cadical_external-prefix/src/cadical_external/src
)

# THIS LINE IS ALSO IMPORTANT
#add_dependencies(talisman flint_external)
add_dependencies(talisman abc_external)


target_include_directories(talisman PRIVATE
  ${CMAKE_BINARY_DIR}/flint_external-prefix/src/flint_external
)



target_compile_options(talisman PRIVATE
    -Wall -Wextra -O3 -DNDEBUG
    -DHAVEGETRUSAGE
    -DHAVEUNLOCKEDIO
)

# ============================================================
# ABC (inside includes/)
# ============================================================

set(ABC_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/includes/abc")

if(NOT EXISTS "${ABC_ROOT}")
  message(FATAL_ERROR "ABC directory not found at includes/abc")
endif()

add_library(abc_stubs OBJECT src/abc_stubs.c)
set_source_files_properties(src/abc_stubs.c PROPERTIES LANGUAGE C)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(abc_stubs PRIVATE -w)
endif()

target_include_directories(talisman PRIVATE
    "${ABC_ROOT}/src"
)

target_compile_definitions(talisman PRIVATE
    USE_ABC=1
    ABC_USE_STDINT_H=1
)

# ============================================================
# Link Everything
# ============================================================

target_link_libraries(talisman PRIVATE
    flint
    cadical
    PBlib::pblib
    cmsgen
    abc
    $<TARGET_OBJECTS:abc_stubs>
    pthread
    dl
    m
)

