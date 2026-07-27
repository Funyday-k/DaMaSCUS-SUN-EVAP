if(NOT DEFINED LIBPHYSICA_SOURCE_DIR)
  message(FATAL_ERROR "LIBPHYSICA_SOURCE_DIR is required")
endif()
if(NOT DEFINED PATCH_SCRIPT)
  message(FATAL_ERROR "PATCH_SCRIPT is required")
endif()
if(NOT DEFINED TEST_OUTPUT_DIR)
  message(FATAL_ERROR "TEST_OUTPUT_DIR is required")
endif()

set(NUMERICS_SOURCE "${LIBPHYSICA_SOURCE_DIR}/src/Numerics.cpp")
set(ORIGINAL_EXPRESSION
  "double inte = prefactor * (a[j] * pow((x - x_j), 3.0) + b[j] * pow((x - x_j), 2.0) + c[j] * (x - x_j) + d[j]);")
set(PATCHED_EXPRESSION
  "double inte = prefactor * (((a[j] * (x - x_j) + b[j]) * (x - x_j) + c[j]) * (x - x_j) + d[j]);")

file(SHA256 "${NUMERICS_SOURCE}" SOURCE_HASH_BEFORE)
foreach(PATCH_ENABLED IN ITEMS ON OFF)
  set(OUTPUT_SOURCE "${TEST_OUTPUT_DIR}/Numerics-${PATCH_ENABLED}.cpp")
  execute_process(
    COMMAND ${CMAKE_COMMAND}
      -DLIBPHYSICA_SOURCE_DIR=${LIBPHYSICA_SOURCE_DIR}
      -DPATCHED_OUTPUT=${OUTPUT_SOURCE}
      -DPATCH_ENABLED=${PATCH_ENABLED}
      -P ${PATCH_SCRIPT}
    RESULT_VARIABLE PATCH_RESULT)
  if(NOT PATCH_RESULT EQUAL 0)
    message(FATAL_ERROR "patch script failed with PATCH_ENABLED=${PATCH_ENABLED}")
  endif()
  file(READ "${OUTPUT_SOURCE}" OUTPUT_CONTENTS)
  if(PATCH_ENABLED)
    string(FIND "${OUTPUT_CONTENTS}" "${PATCHED_EXPRESSION}" EXPECTED_POSITION)
  else()
    string(FIND "${OUTPUT_CONTENTS}" "${ORIGINAL_EXPRESSION}" EXPECTED_POSITION)
  endif()
  if(EXPECTED_POSITION EQUAL -1)
    message(FATAL_ERROR
      "generated source did not match PATCH_ENABLED=${PATCH_ENABLED}")
  endif()
endforeach()

file(SHA256 "${NUMERICS_SOURCE}" SOURCE_HASH_AFTER)
if(NOT SOURCE_HASH_BEFORE STREQUAL SOURCE_HASH_AFTER)
  message(FATAL_ERROR "patch script modified the fetched libphysica source")
endif()
