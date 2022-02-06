cmake_minimum_required(VERSION 3.18)

project(GeodesicSim LANGUAGES CXX)

# ---------------------------------------------
# macOS fixes: make sure <vector>, <string>, etc. are found
if (APPLE)
    message(STATUS "Detected macOS - setting up standard library paths")

    execute_process(
        COMMAND xcrun --show-sdk-path
        OUTPUT_VARIABLE MACOSX_SYSROOT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    set(CMAKE_OSX_SYSROOT ${MACOSX_SYSROOT})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isysroot ${MACOSX_SYSROOT} -stdlib=libc++")
endif()
# ---------------------------------------------

option(USE_CUDA "Enable CUDA support" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Set up threading manually
set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)

# Find dependencies
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)

# Common sources
set(SOURCES_COMMON
    src/file_io.cpp
    src/cpu_sim.cpp
)

# Optional CUDA sources
if (USE_CUDA)
    enable_language(CUDA)
    list(APPEND SOURCES_COMMON src/cuda_sim.cu)
endif()

# Main executable
add_executable(geodesic_gl
    src/main.cpp
    src/simulation_engine.cpp
    src/opengl_utils.cpp
    src/renderer.cpp
    ${SOURCES_COMMON}
)

# Link libraries
target_link_libraries(geodesic_gl
    PRIVATE
        ${CMAKE_THREAD_LIBS_INIT}
        OpenGL::GL
        GLEW::GLEW
        glfw
)

# If CUDA enabled
if (USE_CUDA)
    target_compile_definitions(geodesic_gl PRIVATE USE_CUDA=1)
    target_link_libraries(geodesic_gl PRIVATE CUDA::cudart)
endif()

# Output binary into bin/ folder
set_target_properties(geodesic_gl PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

file(COPY shaders DESTINATION ${CMAKE_BINARY_DIR}/bin)