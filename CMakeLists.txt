cmake_minimum_required(VERSION 3.12)

project(DaMaSCUS-SUN VERSION "0.1.1")

include(GNUInstallDirs)

# Default to Release on single-config generators so local builds use optimization.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Keep older vendored dependencies configurable under newer CMake.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Minimum policy version for vendored deps")
if(POLICY CMP0167)
  # Keep FindBoost support for distributions that do not ship BoostConfig.cmake.
  cmake_policy(SET CMP0167 OLD)
  set(CMAKE_POLICY_DEFAULT_CMP0167 OLD)
endif()
if(POLICY CMP0169)
  # Manual population is needed to apply the small obscura compatibility patch
  # before add_subdirectory(), while retaining compatibility with CMake 3.12.
  cmake_policy(SET CMP0169 OLD)
  set(CMAKE_POLICY_DEFAULT_CMP0169 OLD)
endif()

# Boost.Math 1.89+ requires C++14. C++14 remains compatible with the pinned
# obscura/libphysica releases and avoids a compiler failure with current Boost.
set(CMAKE_CXX_STANDARD 14)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

# The Monte Carlo hot path crosses static-library boundaries on every integrator
# step (this project -> libphysica interpolation and linear algebra), so without
# link-time optimization none of it can be inlined across translation units.
# Runtime dominates build time for production runs, so this defaults to on.
# Set these before the dependencies are added so their targets inherit them.
option(DAMASCUS_ENABLE_LTO "Enable link-time optimization for optimized builds" ON)
option(DAMASCUS_ENABLE_NATIVE_ARCH "Tune the build for the building machine's CPU" OFF)
option(DAMASCUS_PATCH_LIBPHYSICA_INTERPOLATION
  "Build libphysica's cubic interpolation in Horner form" ON)

if(DAMASCUS_ENABLE_LTO AND NOT CODE_COVERAGE)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT DAMASCUS_IPO_SUPPORTED OUTPUT DAMASCUS_IPO_OUTPUT)
  if(DAMASCUS_IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
    message(STATUS "Link-time optimization enabled for optimized builds")
  else()
    message(STATUS "Link-time optimization unavailable: ${DAMASCUS_IPO_OUTPUT}")
  endif()
endif()

# Opt-in only: a native build must run on the machine that compiled it. This is
# safe for the Slurm workflow in bin/, which compiles on the target node.
if(DAMASCUS_ENABLE_NATIVE_ARCH)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-march=native DAMASCUS_HAS_MARCH_NATIVE)
  if(DAMASCUS_HAS_MARCH_NATIVE)
    add_compile_options(-march=native)
  else()
    check_cxx_compiler_flag(-mcpu=native DAMASCUS_HAS_MCPU_NATIVE)
    if(DAMASCUS_HAS_MCPU_NATIVE)
      add_compile_options(-mcpu=native)
    else()
      message(WARNING
        "DAMASCUS_ENABLE_NATIVE_ARCH is ON, but this compiler supports neither "
        "-march=native nor -mcpu=native; building without CPU tuning.")
    endif()
  endif()
endif()

set(GENERATED_DIR ${PROJECT_BINARY_DIR}/generated)
set(INCLUDE_DIR   ${PROJECT_SOURCE_DIR}/include)
set(SRC_DIR       ${PROJECT_SOURCE_DIR}/src)
set(TESTS_DIR     ${PROJECT_SOURCE_DIR}/tests)

# Runtime data is installed relative to the executable so a complete install
# prefix can be moved or copied without rebuilding. Keep the path relative: it
# also makes `cmake --install --prefix <dir>` work as expected.
set(DAMASCUS_INSTALL_DATADIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
if(IS_ABSOLUTE "${CMAKE_INSTALL_BINDIR}" OR IS_ABSOLUTE "${CMAKE_INSTALL_DATADIR}")
  message(FATAL_ERROR
    "CMAKE_INSTALL_BINDIR and CMAKE_INSTALL_DATADIR must be relative so the "
    "installed directory remains relocatable.")
endif()
file(RELATIVE_PATH DAMASCUS_DATA_PATH_FROM_BINDIR
  "/${CMAKE_INSTALL_BINDIR}"
  "/${DAMASCUS_INSTALL_DATADIR}/model_agss09.dat")

# Mirror runtime data into the build tree. Build-tree executables and tests use
# the same executable-relative lookup as installed binaries, without embedding
# an absolute checkout path.
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/${DAMASCUS_INSTALL_DATADIR}")
configure_file(
  "${PROJECT_SOURCE_DIR}/data/model_agss09.dat"
  "${PROJECT_BINARY_DIR}/${DAMASCUS_INSTALL_DATADIR}/model_agss09.dat"
  COPYONLY)

#External projects
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Threads REQUIRED)
find_package(Boost 1.65 REQUIRED)

# obscura's v1.0.1 build uses these cache variable names. Resolve them here so
# missing libconfig++ produces a useful configure-time error on every platform.
find_path(LIBCONFIG_INCLUDE_DIRs NAMES libconfig.h++)
find_library(LIBCONFIGPP_LIBRARY NAMES config++)
if(NOT LIBCONFIG_INCLUDE_DIRs OR NOT LIBCONFIGPP_LIBRARY)
  message(FATAL_ERROR
    "libconfig++ was not found. Install libconfig++ development files "
    "(Ubuntu/Debian: libconfig++-dev; Homebrew: libconfig), or set "
    "LIBCONFIG_INCLUDE_DIRs and LIBCONFIGPP_LIBRARY explicitly.")
endif()

include(FetchContent)
# Declare googletest before adding dependencies. FetchContent keeps the first
# declaration for a content name; obscura/libphysica also declare
# "googletest", so declaring it here ensures the source lives at the root
# project's stable path on clean and repeated configurations.
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  # Exact commit for release-1.10.0. A commit is reproducible even if a tag is
  # moved upstream.
  GIT_TAG        703bd9caab50b139428cea1aaff9974ebee5742e
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      "" )

# Declare libphysica before obscura. FetchContent keeps the first declaration,
# so this pins obscura's transitive dependency as well.
FetchContent_Declare(libphysica
  GIT_REPOSITORY https://github.com/temken/libphysica.git
  # Exact commit for v0.1.2.
  GIT_TAG        fefbe47993bbd343277718708407fa7a40629738
  # obscura v1.0.1 calls add_subdirectory() without an explicit binary dir, so
  # its dependency must remain inside the obscura source tree.
  SOURCE_DIR     "${FETCHCONTENT_BASE_DIR}/obscura-src/external/libphysica")

# obscura
FetchContent_Declare(obscura
  GIT_REPOSITORY  https://github.com/temken/obscura.git
  # Exact commit for v1.0.1.
  GIT_TAG         4b1b9d87f8a69da4d8081175d322b66e48314328
  CMAKE_ARGS     "-DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF" )
FetchContent_GetProperties(obscura)
if(NOT obscura_POPULATED)
    FetchContent_Populate(obscura)
    execute_process(
      COMMAND ${CMAKE_COMMAND}
        -DOBSCURA_SOURCE_DIR=${obscura_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/PatchObscuraCrystalInitialization.cmake
      RESULT_VARIABLE OBSCURA_PATCH_RESULT)
    if(NOT OBSCURA_PATCH_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to patch obscura crystal detector initialization")
    endif()
    # The dependency is linked into our executable, but its own executables,
    # tests, and source-tree install rules are not part of this deployment.
    add_subdirectory(${obscura_SOURCE_DIR} ${obscura_BINARY_DIR} EXCLUDE_FROM_ALL)
    # libphysica is populated by obscura's own FetchContent call inside the
    # add_subdirectory() above, so its source directory is only known here.
    # Generate a build-tree copy of Numerics.cpp and replace that source in the
    # target. This keeps the fetched checkout pristine and makes toggling the
    # optimization in either direction deterministic.
    FetchContent_GetProperties(libphysica)
    get_target_property(LIBPHYSICA_TARGET_SOURCE_DIR libphysica SOURCE_DIR)
    if(NOT libphysica_SOURCE_DIR)
      get_filename_component(libphysica_SOURCE_DIR "${LIBPHYSICA_TARGET_SOURCE_DIR}" DIRECTORY)
    endif()
    set(DAMASCUS_LIBPHYSICA_NUMERICS_BUILD_SOURCE
      "${PROJECT_BINARY_DIR}/patched_dependencies/libphysica/Numerics.cpp")
    execute_process(
      COMMAND ${CMAKE_COMMAND}
        -DLIBPHYSICA_SOURCE_DIR=${libphysica_SOURCE_DIR}
        -DPATCHED_OUTPUT=${DAMASCUS_LIBPHYSICA_NUMERICS_BUILD_SOURCE}
        -DPATCH_ENABLED=${DAMASCUS_PATCH_LIBPHYSICA_INTERPOLATION}
        -P ${PROJECT_SOURCE_DIR}/cmake/PatchLibphysicaInterpolation.cmake
      RESULT_VARIABLE LIBPHYSICA_PATCH_RESULT)
    if(NOT LIBPHYSICA_PATCH_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to prepare libphysica cubic interpolation source")
    endif()

    get_target_property(LIBPHYSICA_SOURCES libphysica SOURCES)
    set(DAMASCUS_LIBPHYSICA_SOURCES)
    set(DAMASCUS_REPLACED_LIBPHYSICA_NUMERICS OFF)
    foreach(LIBPHYSICA_SOURCE IN LISTS LIBPHYSICA_SOURCES)
      if(IS_ABSOLUTE "${LIBPHYSICA_SOURCE}")
        set(LIBPHYSICA_ABSOLUTE_SOURCE "${LIBPHYSICA_SOURCE}")
      else()
        set(LIBPHYSICA_ABSOLUTE_SOURCE
          "${LIBPHYSICA_TARGET_SOURCE_DIR}/${LIBPHYSICA_SOURCE}")
      endif()
      get_filename_component(LIBPHYSICA_SOURCE_NAME
        "${LIBPHYSICA_ABSOLUTE_SOURCE}" NAME)
      if(LIBPHYSICA_SOURCE_NAME STREQUAL "Numerics.cpp")
        list(APPEND DAMASCUS_LIBPHYSICA_SOURCES
          "${DAMASCUS_LIBPHYSICA_NUMERICS_BUILD_SOURCE}")
        set(DAMASCUS_REPLACED_LIBPHYSICA_NUMERICS ON)
      else()
        list(APPEND DAMASCUS_LIBPHYSICA_SOURCES
          "${LIBPHYSICA_ABSOLUTE_SOURCE}")
      endif()
    endforeach()
    if(NOT DAMASCUS_REPLACED_LIBPHYSICA_NUMERICS)
      message(FATAL_ERROR "libphysica target no longer contains Numerics.cpp")
    endif()
    set_property(TARGET libphysica PROPERTY SOURCES
      "${DAMASCUS_LIBPHYSICA_SOURCES}")
    # Both pinned dependencies declare C++11 internally. Override their target
    # properties because current Boost.Math headers require C++14.
    set_target_properties(libphysica libobscura PROPERTIES
      CXX_STANDARD 14
      CXX_STANDARD_REQUIRED ON
      CXX_EXTENSIONS OFF)
endif()

# Refresh the build-tree version header on every build.
set(VERSION_HEADER_TEMPLATE ${INCLUDE_DIR}/version.hpp.in)
set(VERSION_HEADER ${GENERATED_DIR}/version.hpp)
add_custom_target(generate_version_header ALL
  COMMAND ${CMAKE_COMMAND}
    "-DSOURCE_DIR=${PROJECT_SOURCE_DIR}"
    "-DTEMPLATE=${VERSION_HEADER_TEMPLATE}"
    "-DOUTPUT=${VERSION_HEADER}"
    "-DPROJECT_NAME=${PROJECT_NAME}"
    "-DPROJECT_VERSION=${PROJECT_VERSION}"
    "-DPROJECT_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}"
    "-DPROJECT_VERSION_MINOR=${PROJECT_VERSION_MINOR}"
    "-DPROJECT_VERSION_PATCH=${PROJECT_VERSION_PATCH}"
    "-DINSTALL_DATA_PATH_FROM_BINDIR=${DAMASCUS_DATA_PATH_FROM_BINDIR}"
    -P ${PROJECT_SOURCE_DIR}/cmake/GenerateVersionHeader.cmake
  BYPRODUCTS ${VERSION_HEADER}
  COMMENT "Refreshing generated version header"
  VERBATIM)

# Source and include directories
include_directories(${INCLUDE_DIR})
add_subdirectory( ${SRC_DIR} )

install(FILES ${PROJECT_SOURCE_DIR}/data/model_agss09.dat
  DESTINATION ${DAMASCUS_INSTALL_DATADIR})
install(FILES ${PROJECT_SOURCE_DIR}/examples/quickstart.cfg
  DESTINATION ${DAMASCUS_INSTALL_DATADIR}/examples)

# Code Coverage Configuration
if(NOT TARGET coverage_config)
  add_library(coverage_config INTERFACE)
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE)
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE)

option (BUILD_TESTING "Build the testing tree." ON)
# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
  enable_testing()
  FetchContent_GetProperties(googletest)
  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(
      ${googletest_SOURCE_DIR} 
      ${googletest_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif()
  add_subdirectory(${TESTS_DIR})
endif()
