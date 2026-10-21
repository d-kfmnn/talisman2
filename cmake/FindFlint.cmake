include(FindPackageHandleStandardArgs)

# Try to find libraries
find_library(FLINT_C_LIBRARIES
  NAMES flint
  DOC "FLINT C libraries"
)

# Try to find headers
find_path(FLINT_C_INCLUDES
  NAMES flint.h
  PATH_SUFFIXES flint
  DOC "Flint C header"
)

find_package_handle_standard_args(Flint
	REQUIRED_VARS FLINT_C_LIBRARIES FLINT_C_INCLUDES)

if (FLINT_FOUND)
  set(FLINT_INCLUDE_DIRS "${FLINT_C_INCLUDES}")

  if (NOT TARGET Flint::Flint)
    add_library(Flint::Flint UNKNOWN IMPORTED)
    set_target_properties(Flint::Flint PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${FLINT_C_INCLUDES}"
      IMPORTED_LOCATION "${FLINT_C_LIBRARIES}")
  endif()
endif()
