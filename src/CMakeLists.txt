# Executable
add_executable(DaMaSCUS-SUN
	main.cpp)

target_compile_options(DaMaSCUS-SUN PUBLIC -Wall -pedantic)

target_link_libraries(DaMaSCUS-SUN
	PUBLIC
		lib_damascus_sun)

target_include_directories(DaMaSCUS-SUN
	PRIVATE
		${GENERATED_DIR} )

install(TARGETS DaMaSCUS-SUN
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(replay-invalid-trajectory
	replay_invalid_trajectory.cpp)

target_compile_options(replay-invalid-trajectory PUBLIC -Wall -pedantic)

target_link_libraries(replay-invalid-trajectory
	PUBLIC
		lib_damascus_sun)

target_include_directories(replay-invalid-trajectory
	PRIVATE
		${GENERATED_DIR})

install(TARGETS replay-invalid-trajectory
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Static library
add_library(lib_damascus_sun STATIC
	Dark_Photon.cpp
	Data_Generation.cpp
	Parameter_Scan.cpp
	Reflection_Spectrum.cpp
	Snapshot_Heartbeat.cpp
	Snapshot_IO.cpp
	Snapshot_Shared_State.cpp
	Simulation_Trajectory.cpp
	Simulation_Utilities.cpp
	Earth_Model.cpp
	Solar_Model.cpp )

set_target_properties(lib_damascus_sun PROPERTIES PREFIX "")

target_compile_options(lib_damascus_sun PUBLIC -Wall -pedantic)

target_include_directories(lib_damascus_sun 
    PRIVATE
		${GENERATED_DIR})

add_dependencies(DaMaSCUS-SUN generate_version_header)
add_dependencies(replay-invalid-trajectory generate_version_header)
add_dependencies(lib_damascus_sun generate_version_header)

target_link_libraries(lib_damascus_sun
	    PUBLIC
		coverage_config
			libobscura
			MPI::MPI_CXX
			Threads::Threads )
