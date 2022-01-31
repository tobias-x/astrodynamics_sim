cmake_minimum_required(VERSION 3.18)

# Project name and version
project(GeodesicSim LANGUAGES CXX)

# Option to enable CUDA support
option(USE_CUDA "Enable CUDA support" OFF)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Common source files
set(SOURCES_COMMON
    src/file_io.cpp
    src/run_timestep_sim.cpp
    src/run_realtime_sim.cpp
    src/cpu_sim.cpp
)

# Handle CUDA support
if (USE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    list(APPEND SOURCES_COMMON src/cuda_sim.cu)
endif()

# Timestep executable
if (USE_CUDA)
    add_executable(timestep_sim_gpu main/timestep_main.cpp ${SOURCES_COMMON})
    target_compile_definitions(timestep_sim_gpu PRIVATE USE_CUDA=1)
    target_link_libraries(timestep_sim_gpu CUDA::cudart)
    set_target_properties(timestep_sim_gpu PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
else()
    add_executable(timestep_sim main/timestep_main.cpp ${SOURCES_COMMON})
    target_compile_definitions(timestep_sim PRIVATE USE_CUDA=0)
    set_target_properties(timestep_sim PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endif()

# Runtime executable
if (USE_CUDA)
    add_executable(runtime_sim_gpu main/runtime_main.cpp ${SOURCES_COMMON})
    target_compile_definitions(runtime_sim_gpu PRIVATE USE_CUDA=1)
    target_link_libraries(runtime_sim_gpu CUDA::cudart)
    set_target_properties(runtime_sim_gpu PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
else()
    add_executable(runtime_sim main/runtime_main.cpp ${SOURCES_COMMON})
    target_compile_definitions(runtime_sim PRIVATE USE_CUDA=0)
    set_target_properties(runtime_sim PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endif()

# Print a summary
if (USE_CUDA)
    message(STATUS "Building with CUDA support: ON")
else()
    message(STATUS "Building with CUDA support: OFF")
endif()
