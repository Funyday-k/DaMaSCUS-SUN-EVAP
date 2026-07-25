# 1. Test executable
add_test(NAME executable_test
	COMMAND DaMaSCUS-SUN ${PROJECT_SOURCE_DIR}/tests/config_executable_smoke.cfg
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME test_Invalid_Snapshot_Interval_Config
	COMMAND ${CMAKE_COMMAND}
		-DPROGRAM=$<TARGET_FILE:DaMaSCUS-SUN>
		-DCONFIG=${PROJECT_SOURCE_DIR}/tests/config_snapshot_invalid_fractional.cfg
		"-DEXPECTED_TEXT='snapshot_interval' must be a positive integer number of seconds."
		-P ${PROJECT_SOURCE_DIR}/cmake/ExpectCommandFailure.cmake
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(test_Invalid_Snapshot_Interval_Config PROPERTIES
	TIMEOUT 8
	LABELS "snapshot;config")

# 2. Unit tests: Dark_Photon.xpp	
add_executable(test_Dark_Photon
	test_Dark_Photon.cpp)
target_include_directories(test_Dark_Photon PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Dark_Photon
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Dark_Photon PUBLIC -Wall -pedantic)
add_test(NAME test_Dark_Photon COMMAND test_Dark_Photon
	WORKING_DIRECTORY ${TESTS_DIR})

# 3. Unit tests: Data_Generation.xpp	
add_executable(test_Data_Generation
	test_Data_Generation.cpp)
target_include_directories(test_Data_Generation PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Data_Generation
	PRIVATE
		lib_damascus_sun
		gtest)
target_compile_options(test_Data_Generation PUBLIC -Wall -pedantic)
add_test(NAME test_Data_Generation COMMAND test_Data_Generation
	WORKING_DIRECTORY ${TESTS_DIR})

# 4. Unit tests: Parameter_Scan.xpp	
add_executable(test_Parameter_Scan
	test_Parameter_Scan.cpp)
target_include_directories(test_Parameter_Scan PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Parameter_Scan
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Parameter_Scan PUBLIC -Wall -pedantic)
add_test(NAME test_Parameter_Scan COMMAND test_Parameter_Scan
	WORKING_DIRECTORY ${TESTS_DIR})

# 5. Unit tests: Reflection_Spectrum.xpp	
add_executable(test_Reflection_Spectrum
	test_Reflection_Spectrum.cpp)
target_include_directories(test_Reflection_Spectrum PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Reflection_Spectrum
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Reflection_Spectrum PUBLIC -Wall -pedantic)
add_test(NAME test_Reflection_Spectrum COMMAND test_Reflection_Spectrum
	WORKING_DIRECTORY ${TESTS_DIR})

# 6. Unit tests: Simulation_Trajectory.xpp	
add_executable(test_Simulation_Trajectory
	test_Simulation_Trajectory.cpp)
target_include_directories(test_Simulation_Trajectory PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Simulation_Trajectory
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Simulation_Trajectory PUBLIC -Wall -pedantic)
add_test(NAME test_Simulation_Trajectory COMMAND test_Simulation_Trajectory
	WORKING_DIRECTORY ${TESTS_DIR})

# 7. Unit tests: Simulation_Utilties.xpp	
add_executable(test_Simulation_Utilities
	test_Simulation_Utilities.cpp)
target_include_directories(test_Simulation_Utilities PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Simulation_Utilities
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Simulation_Utilities PUBLIC -Wall -pedantic)
add_test(NAME test_Simulation_Utilities COMMAND test_Simulation_Utilities
	WORKING_DIRECTORY ${TESTS_DIR})

# 7b. Physics-facing validation: analytic limits, invariants, and sampling laws.
add_executable(test_Physics_Validation
	test_Physics_Validation.cpp)
target_include_directories(test_Physics_Validation PRIVATE ${GENERATED_DIR})
target_link_libraries(test_Physics_Validation
	PRIVATE
		lib_damascus_sun
		gtest_main)
target_compile_options(test_Physics_Validation PUBLIC -Wall -pedantic)
add_test(NAME test_Physics_Validation COMMAND test_Physics_Validation
	WORKING_DIRECTORY ${TESTS_DIR})
set_tests_properties(test_Physics_Validation PROPERTIES
	TIMEOUT 30
	LABELS "physics-validation;unit")

# 8. Unit tests: Solar_Model.xpp	
add_executable(test_Solar_Model
	test_Solar_Model.cpp)
target_include_directories(test_Solar_Model PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Solar_Model
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Solar_Model PUBLIC -Wall -pedantic)
add_test(NAME test_Solar_Model COMMAND test_Solar_Model
	WORKING_DIRECTORY ${TESTS_DIR})
if(MPIEXEC_EXECUTABLE
	AND MPIEXEC_NUMPROC_FLAG
	AND DEFINED MPIEXEC_MAX_NUMPROCS
	AND MPIEXEC_MAX_NUMPROCS GREATER_EQUAL 3)
	add_test(NAME test_Solar_Model_MPI_3Ranks
		COMMAND
			${MPIEXEC_EXECUTABLE}
			${MPIEXEC_NUMPROC_FLAG} 3
			${MPIEXEC_PREFLAGS}
			$<TARGET_FILE:test_Solar_Model>
			--gtest_filter=TestSolarModel.TestTotalDMScatteringRateRegularGridNodes
			${MPIEXEC_POSTFLAGS}
		WORKING_DIRECTORY ${TESTS_DIR})
	set_tests_properties(test_Solar_Model_MPI_3Ranks PROPERTIES
		PROCESSORS 3
		RUN_SERIAL TRUE
		TIMEOUT 60
		LABELS "solar-model;mpi")
else()
	message(STATUS "Skipping three-rank solar-model grid test: MPI launcher or explicit three-rank capacity is unavailable (MPIEXEC_EXECUTABLE='${MPIEXEC_EXECUTABLE}', MPIEXEC_MAX_NUMPROCS='${MPIEXEC_MAX_NUMPROCS}')")
endif()

# 9. Snapshot serialization and merge tests (no solar model/runtime data)
add_executable(test_Snapshot_IO
	test_Snapshot_IO.cpp)
target_include_directories(test_Snapshot_IO PRIVATE ${GENERATED_DIR})
target_link_libraries(test_Snapshot_IO
	PRIVATE
		lib_damascus_sun
		gtest_main)
target_compile_options(test_Snapshot_IO PUBLIC -Wall -pedantic)
add_test(NAME test_Snapshot_IO COMMAND test_Snapshot_IO
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(test_Snapshot_IO PROPERTIES
	TIMEOUT 10
	LABELS "snapshot;unit")

# 10. Real four-rank heartbeat test with barrier-synchronized local epochs.
add_executable(test_Snapshot_Heartbeat_MPI
	test_Snapshot_Heartbeat_MPI.cpp)
target_include_directories(test_Snapshot_Heartbeat_MPI PRIVATE ${GENERATED_DIR})
target_link_libraries(test_Snapshot_Heartbeat_MPI
	PRIVATE
		lib_damascus_sun
		MPI::MPI_CXX)
target_compile_options(test_Snapshot_Heartbeat_MPI PUBLIC -Wall -pedantic)
if(MPIEXEC_EXECUTABLE
	AND MPIEXEC_NUMPROC_FLAG
	AND DEFINED MPIEXEC_MAX_NUMPROCS
	AND MPIEXEC_MAX_NUMPROCS GREATER_EQUAL 4)
	add_test(NAME test_Snapshot_Heartbeat_MPI_4Ranks
		COMMAND
			${MPIEXEC_EXECUTABLE}
			${MPIEXEC_NUMPROC_FLAG} 4
			${MPIEXEC_PREFLAGS}
			$<TARGET_FILE:test_Snapshot_Heartbeat_MPI>
			${MPIEXEC_POSTFLAGS}
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
	set_tests_properties(test_Snapshot_Heartbeat_MPI_4Ranks PROPERTIES
		PROCESSORS 4
		RUN_SERIAL TRUE
		TIMEOUT 35
		LABELS "snapshot;mpi;slow")
else()
	message(STATUS "Skipping four-rank snapshot heartbeat test: MPI launcher or explicit four-rank capacity is unavailable (MPIEXEC_EXECUTABLE='${MPIEXEC_EXECUTABLE}', MPIEXEC_MAX_NUMPROCS='${MPIEXEC_MAX_NUMPROCS}')")
endif()

# 10b. MPI straggler regression: early ranks keep doing useful tail work.
add_executable(test_MPI_Tail_Work
	test_MPI_Tail_Work.cpp)
target_include_directories(test_MPI_Tail_Work PRIVATE ${GENERATED_DIR})
target_link_libraries(test_MPI_Tail_Work
	PRIVATE
		lib_damascus_sun
		MPI::MPI_CXX)
target_compile_options(test_MPI_Tail_Work PUBLIC -Wall -pedantic)
if(MPIEXEC_EXECUTABLE
	AND MPIEXEC_NUMPROC_FLAG
	AND DEFINED MPIEXEC_MAX_NUMPROCS
	AND MPIEXEC_MAX_NUMPROCS GREATER_EQUAL 2)
	add_test(NAME test_MPI_Tail_Work_2Ranks
		COMMAND
			${MPIEXEC_EXECUTABLE}
			${MPIEXEC_NUMPROC_FLAG} 2
			${MPIEXEC_PREFLAGS}
			$<TARGET_FILE:test_MPI_Tail_Work>
			${MPIEXEC_POSTFLAGS}
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
	set_tests_properties(test_MPI_Tail_Work_2Ranks PROPERTIES
		PROCESSORS 2
		RUN_SERIAL TRUE
		TIMEOUT 10
		LABELS "mpi;performance")
else()
	message(STATUS "Skipping two-rank MPI tail-work test: MPI launcher or explicit two-rank capacity is unavailable (MPIEXEC_EXECUTABLE='${MPIEXEC_EXECUTABLE}', MPIEXEC_MAX_NUMPROCS='${MPIEXEC_MAX_NUMPROCS}')")
endif()

# 11. Python workflow tests.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
	add_test(NAME test_Numerical_Stability_Python
		COMMAND ${Python3_EXECUTABLE} ${TESTS_DIR}/test_numerical_stability.py
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
	set_tests_properties(test_Numerical_Stability_Python PROPERTIES
		TIMEOUT 30
		LABELS "numerical-stability;unit")
	add_test(NAME test_Physics_Validation_Workflow
		COMMAND ${Python3_EXECUTABLE} ${TESTS_DIR}/test_physics_validation.py
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
	set_tests_properties(test_Physics_Validation_Workflow PROPERTIES
		TIMEOUT 10
		LABELS "physics-validation;unit")
else()
	message(STATUS "Skipping Python workflow tests: Python3 interpreter is unavailable")
endif()
