cmake_minimum_required(VERSION 3.16)
project(talisman C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Allow user to specify Cadical include dir if it's not in standard locations
set(CADICAL_INCLUDE_DIR "" CACHE PATH "Path to CaDiCaL include directory")

# Sanitizer option
option(SANITIZE "Enable AddressSanitizer and UBSan" OFF)

# Compile flags
set(COMMON_FLAGS -Wall -Wextra -O3 -DNDEBUG -DHAVEGETRUSAGE -DHAVEUNLOCKEDIO)
if(SANITIZE)
  list(APPEND COMMON_FLAGS -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer -g)
endif()
add_compile_options(${COMMON_FLAGS})

# Source files
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS src/*.cpp)

# Remove testers from main target
list(REMOVE_ITEM SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/cadical_tester.cpp")

# Compile aiger.c as an object
add_library(aiger_obj OBJECT includes/aiger/aiger.c)
set_source_files_properties(includes/aiger/aiger.c PROPERTIES LANGUAGE C)

# Suppress warnings just for aiger.c
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(aiger_obj PRIVATE -w)
elseif(MSVC)
  target_compile_options(aiger_obj PRIVATE /w)
endif()

# Executable
add_executable(talisman ${SRC_FILES} $<TARGET_OBJECTS:aiger_obj>)

# ---- Find GMP and GMPXX manually ----
find_library(GMP_LIBRARY NAMES gmp PATHS ${CMAKE_LIBRARY_PATH})
find_path(GMP_INCLUDE_DIR NAMES gmp.h PATHS ${CMAKE_INCLUDE_PATH})
if(NOT GMP_LIBRARY OR NOT GMP_INCLUDE_DIR)
  message(FATAL_ERROR "GMP library or headers not found")
endif()
message(STATUS "GMP library found at: ${GMP_LIBRARY}")
message(STATUS "GMP include dir found at: ${GMP_INCLUDE_DIR}")

find_library(GMPXX_LIBRARY NAMES gmpxx PATHS ${CMAKE_LIBRARY_PATH})
find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h PATHS ${CMAKE_INCLUDE_PATH})
if(NOT GMPXX_LIBRARY OR NOT GMPXX_INCLUDE_DIR)
  message(FATAL_ERROR "GMPXX library not found")
endif()
message(STATUS "GMPXX library found at: ${GMPXX_LIBRARY}")

# ---- FLINT: force user-provided paths (no find_package) ----
set(FLINT_INCLUDE_DIR "" CACHE PATH "Path to FLINT include directory (parent include dir, contains flint/)")
set(FLINT_LIBRARY "" CACHE FILEPATH "Path to FLINT library file (libflint.so or libflint.a)")

if(NOT FLINT_INCLUDE_DIR OR NOT FLINT_LIBRARY)
  message(FATAL_ERROR "Provide -DFLINT_INCLUDE_DIR=/path/to/include and -DFLINT_LIBRARY=/path/to/libflint.so (or .a)")
endif()

# Validate FLINT include dir: must contain flint/nmod.h
if(NOT EXISTS "${FLINT_INCLUDE_DIR}/flint/nmod.h")
  message(FATAL_ERROR
    "flint/nmod.h not found under ${FLINT_INCLUDE_DIR}. "
    "Set FLINT_INCLUDE_DIR to the parent include dir (e.g. ~/.local/include)."
  )
endif()

# Support both include styles:
#   #include <flint/nmod.h>  -> needs -I${FLINT_INCLUDE_DIR}
#   #include <nmod_mat.h>    -> needs -I${FLINT_INCLUDE_DIR}/flint
set(FLINT_INCLUDE_DIRS
  "${FLINT_INCLUDE_DIR}"
  "${FLINT_INCLUDE_DIR}/flint"
)

# Validate FLINT library file
if(NOT EXISTS "${FLINT_LIBRARY}")
  message(FATAL_ERROR "FLINT library file does not exist: ${FLINT_LIBRARY}")
endif()

# Create an imported library target Flint::Flint
add_library(Flint::Flint UNKNOWN IMPORTED)
set_target_properties(Flint::Flint PROPERTIES
  IMPORTED_LOCATION "${FLINT_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${FLINT_INCLUDE_DIRS}"
)

message(STATUS "Forcing user-provided FLINT")
message(STATUS "  FLINT_INCLUDE_DIR  = ${FLINT_INCLUDE_DIR}")
message(STATUS "  FLINT_INCLUDE_DIRS = ${FLINT_INCLUDE_DIRS}")
message(STATUS "  FLINT_LIBRARY      = ${FLINT_LIBRARY}")

# ---- Find CaDiCaL library and headers ----
set(CADICAL_LIBRARY_DIR "" CACHE PATH "Path to CaDiCaL library directory")

find_library(CADICAL_LIBRARY
  NAMES cadical
  HINTS ${CADICAL_LIBRARY_DIR}
  PATHS ${CMAKE_LIBRARY_PATH}
  NO_DEFAULT_PATH
)

if(NOT CADICAL_LIBRARY)
  message(FATAL_ERROR "CaDiCaL library not found. Please set CADICAL_LIBRARY_DIR by adding \"-DCADICAL_LIBRARY_DIR=/path/to/dir/\"")
endif()

message(STATUS "Using CaDiCaL headers from: ${CADICAL_INCLUDE_DIR}")

if(NOT CADICAL_INCLUDE_DIR)
  message(FATAL_ERROR "cadical headers not found (cadical.hpp)")
endif()

message(STATUS "cadical library found at: ${CADICAL_LIBRARY}")
message(STATUS "cadical include dir found at: ${CADICAL_INCLUDE_DIR}")

# Candidate (half-/full-adder, K-cover) detection is implemented in-tree in
# src/candidates.cpp -- no external EDA library dependency needed.
target_link_libraries(${PROJECT_NAME} PRIVATE
  pthread
  m
)

# ---- Include directories ----
# GMP includes
target_include_directories(talisman PRIVATE ${GMP_INCLUDE_DIR} ${GMPXX_INCLUDE_DIR})
# FLINT includes (support both include styles)
target_include_directories(talisman PRIVATE ${FLINT_INCLUDE_DIRS})
# CaDiCaL includes
target_include_directories(talisman PRIVATE ${CADICAL_INCLUDE_DIR})

# Add pblib subdirectory
add_subdirectory(includes/pblib)
# Suppress warnings from PBLib
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(pblib PRIVATE -w)
elseif(MSVC)
  target_compile_options(pblib PRIVATE /w)
endif()
target_include_directories(talisman PRIVATE includes/pblib)

# cmsgen includes
find_package(cmsgen REQUIRED)
target_include_directories(talisman PRIVATE ${CMSGEN_INCLUDE_DIRS})
target_link_libraries(talisman PRIVATE cmsgen)

# ---- Link Libraries ----
target_link_libraries(talisman
  PRIVATE
    Flint::Flint
    ${GMPXX_LIBRARY}
    ${GMP_LIBRARY}
    ${CADICAL_LIBRARY}
    PBlib::pblib
)

# Apply sanitizer link flags
if(SANITIZE)
  target_link_options(talisman PRIVATE -fsanitize=address -fsanitize=undefined)
endif()
