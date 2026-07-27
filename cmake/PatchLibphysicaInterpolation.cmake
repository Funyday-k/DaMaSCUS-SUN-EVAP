if(NOT DEFINED LIBPHYSICA_SOURCE_DIR)
  message(FATAL_ERROR "LIBPHYSICA_SOURCE_DIR is required")
endif()
if(NOT DEFINED PATCHED_OUTPUT)
  message(FATAL_ERROR "PATCHED_OUTPUT is required")
endif()
if(NOT DEFINED PATCH_ENABLED)
  message(FATAL_ERROR "PATCH_ENABLED is required")
endif()

set(NUMERICS_SOURCE "${LIBPHYSICA_SOURCE_DIR}/src/Numerics.cpp")
if(NOT EXISTS "${NUMERICS_SOURCE}")
  message(FATAL_ERROR "libphysica numerics source not found: ${NUMERICS_SOURCE}")
endif()

# Interpolation::Interpolate is the hottest lookup in the trajectory integrator:
# Solar_Model::Mass() alone is evaluated six times per Runge-Kutta attempt. At
# -O3 without -ffast-math the compiler folds pow(dx, 2.0) into a multiplication
# but keeps pow(dx, 3.0) as a libm call, so the cubic term costs a function call
# per evaluation. Horner's form is the same polynomial without any pow().
set(ORIGINAL_EXPRESSION
  "double inte = prefactor * (a[j] * pow((x - x_j), 3.0) + b[j] * pow((x - x_j), 2.0) + c[j] * (x - x_j) + d[j]);")
set(PATCHED_EXPRESSION
  "double inte = prefactor * (((a[j] * (x - x_j) + b[j]) * (x - x_j) + c[j]) * (x - x_j) + d[j]);")

file(READ "${NUMERICS_SOURCE}" CONTENTS)
string(FIND "${CONTENTS}" "${ORIGINAL_EXPRESSION}" ORIGINAL_POSITION)
string(FIND "${CONTENTS}" "${PATCHED_EXPRESSION}" PATCHED_POSITION)
if(ORIGINAL_POSITION EQUAL -1 AND PATCHED_POSITION EQUAL -1)
  message(FATAL_ERROR
    "unsupported libphysica Numerics.cpp: Interpolation::Interpolate body changed")
endif()

if(PATCH_ENABLED)
  string(REPLACE "${ORIGINAL_EXPRESSION}" "${PATCHED_EXPRESSION}" OUTPUT_CONTENTS "${CONTENTS}")
  message(STATUS "Using Horner-form libphysica cubic interpolation")
else()
  string(REPLACE "${PATCHED_EXPRESSION}" "${ORIGINAL_EXPRESSION}" OUTPUT_CONTENTS "${CONTENTS}")
  message(STATUS "Using upstream libphysica cubic interpolation")
endif()

get_filename_component(PATCHED_OUTPUT_DIR "${PATCHED_OUTPUT}" DIRECTORY)
file(MAKE_DIRECTORY "${PATCHED_OUTPUT_DIR}")
if(EXISTS "${PATCHED_OUTPUT}")
  file(READ "${PATCHED_OUTPUT}" EXISTING_OUTPUT_CONTENTS)
  if(EXISTING_OUTPUT_CONTENTS STREQUAL OUTPUT_CONTENTS)
    return()
  endif()
endif()
file(WRITE "${PATCHED_OUTPUT}" "${OUTPUT_CONTENTS}")
