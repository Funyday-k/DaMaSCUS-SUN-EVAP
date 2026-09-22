# Isolate each smoke-test invocation from the previous run's retained results.
if(NOT DEFINED PROGRAM OR NOT DEFINED CONFIG OR NOT DEFINED TEST_OUTPUT_DIR)
  message(FATAL_ERROR "PROGRAM, CONFIG and TEST_OUTPUT_DIR are required")
endif()
file(REMOVE_RECURSE "${TEST_OUTPUT_DIR}")
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")
execute_process(COMMAND "${PROGRAM}" "${CONFIG}"
  WORKING_DIRECTORY "${TEST_OUTPUT_DIR}" TIMEOUT 20
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Quickstart failed: ${result}\n${output}\n${error}")
endif()
