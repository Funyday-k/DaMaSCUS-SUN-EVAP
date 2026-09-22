if(NOT DEFINED BUILD_DIR OR NOT DEFINED TEST_OUTPUT_DIR)
  message(FATAL_ERROR "BUILD_DIR and TEST_OUTPUT_DIR are required")
endif()
set(prefix "${TEST_OUTPUT_DIR}/prefix")
set(relocated "${TEST_OUTPUT_DIR}/relocated")
file(REMOVE_RECURSE "${TEST_OUTPUT_DIR}")
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}/working")
execute_process(COMMAND "${CMAKE_COMMAND}" --install "${BUILD_DIR}" --config "${BUILD_CONFIG}" --prefix "${prefix}"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Installation failed: ${output}\n${error}")
endif()
file(RENAME "${prefix}" "${relocated}")
if(NOT EXISTS "${relocated}/${DATA_DIR}/model_agss09.dat")
  message(FATAL_ERROR "Installed solar model is missing")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E env
  --unset=DAMASCUS_SUN_SOLAR_MODEL --unset=DAMASCUS_SUN_DATA_DIR
  "${relocated}/${BIN_DIR}/DaMaSCUS-SUN"
  "${relocated}/${DATA_DIR}/examples/quickstart.cfg"
  WORKING_DIRECTORY "${TEST_OUTPUT_DIR}/working" TIMEOUT 15
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0 OR NOT output MATCHES "Finished in")
  message(FATAL_ERROR "Relocated quickstart failed: ${result}\n${output}\n${error}")
endif()
if(NOT EXISTS "${TEST_OUTPUT_DIR}/working/quickstart_results/results_-2.000000_-32.000000/bincount.tsv")
  message(FATAL_ERROR "Relocated quickstart did not produce bincount.tsv")
endif()
