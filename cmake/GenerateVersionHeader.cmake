# Regenerate the build-tree version header so rebuilt binaries report the
# current branch and commit.
if(NOT DEFINED SOURCE_DIR OR NOT DEFINED TEMPLATE OR NOT DEFINED OUTPUT
   OR NOT DEFINED INSTALL_DATA_PATH_FROM_BINDIR)
  message(FATAL_ERROR
    "SOURCE_DIR, TEMPLATE, OUTPUT, and INSTALL_DATA_PATH_FROM_BINDIR must be defined")
endif()

execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY "${SOURCE_DIR}"
  RESULT_VARIABLE git_branch_result
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT git_branch_result EQUAL 0 OR GIT_BRANCH STREQUAL "")
  set(GIT_BRANCH "unknown")
endif()

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY "${SOURCE_DIR}"
  RESULT_VARIABLE git_commit_result
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT git_commit_result EQUAL 0 OR GIT_COMMIT_HASH STREQUAL "")
  set(GIT_COMMIT_HASH "unknown")
endif()

execute_process(
  COMMAND git status --porcelain --untracked-files=normal
  WORKING_DIRECTORY "${SOURCE_DIR}"
  RESULT_VARIABLE git_status_result
  OUTPUT_VARIABLE git_status_output
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(git_status_result EQUAL 0 AND NOT git_status_output STREQUAL "")
  string(APPEND GIT_COMMIT_HASH "-dirty")
endif()

set(PROJECT_SOURCE_DIR "${SOURCE_DIR}")
get_filename_component(output_dir "${OUTPUT}" DIRECTORY)
file(MAKE_DIRECTORY "${output_dir}")
configure_file("${TEMPLATE}" "${OUTPUT}" @ONLY)
