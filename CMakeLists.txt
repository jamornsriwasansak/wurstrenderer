cmake_minimum_required(VERSION 2.8)
project(wurstrenderer)

# C++ Mode
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -DAVX2")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

# output directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sandbox/archive")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sandbox/library")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sandbox")

# embree
set(EMBREE_ISA_AVX OFF CACHE BOOL "Enables AVX ISA" FORCE)
set(EMBREE_ISA_AVX2 ON CACHE BOOL "Enables AVX2 ISA" FORCE)
set(EMBREE_ISA_SSE2 OFF CACHE BOOL "Enables SSE ISA" FORCE)
set(EMBREE_ISA_SSE42 OFF CACHE BOOL "Enables SSE4.2 ISA" FORCE)
set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "Build Embree with support for ISPC applications." FORCE)
set(EMBREE_GEOMETRY_CURVE OFF CACHE BOOL "Enables support for curve geometries" FORCE)
set(EMBREE_GEOMETRY_GRID OFF CACHE BOOL "Enables support for grid geometries" FORCE)
set(EMBREE_GEOMETRY_INSTANCE OFF CACHE BOOL "Enables support for instances" FORCE)
set(EMBREE_GEOMETRY_POINT OFF CACHE BOOL "Enables support for point geometries" FORCE)
set(EMBREE_GEOMETRY_QUAD OFF CACHE BOOL "Enables support for quad geometries" FORCE)
set(EMBREE_GEOMETRY_SUBDIVISION OFF CACHE BOOL "Enables support for subdiv geometries" FORCE)
set(EMBREE_RAY_PACKETS OFF CACHE BOOL "Enabled support for ray packets" FORCE)
set(EMBREE_STAT_COUNTERS OFF CACHE BOOL "Enables statistic counters." FORCE)
set(EMBREE_STATIC_LIB ON CACHE BOOL "Build Embree as a static library." FORCE)
set(EMBREE_TASKING_SYSTEM "INTERNAL" CACHE STRING "Selects tasking system" FORCE)
set(EMBREE_TESTING_INTENSITY "0" CACHE BOOL "Intensity of testing (0 = no testing, 1 = verify and tutorials, 2 = light testing, 3 = intensive testing" FORCE)
set(EMBREE_TUTORIALS OFF CACHE BOOL "Enable to build Embree tutorials" FORCE)
add_subdirectory(dependencies/embree)
include_directories(dependencies/embree/include)

# assimp
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF CACHE INTERNAL "If the supplementary tools for Assimp are built in addition to the library.")
set(ASSIMP_BUILD_TESTS OFF CACHE INTERNAL "If the test suite for Assimp is built in addition to the library.")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build package with shared libraries.")
set(ASSIMP_BUILD_STATIC_LIB ON CACHE BOOL "Build a static version of the library")
add_subdirectory(dependencies/assimp)
include_directories(dependencies/assimp/include)

# pbrt-parser
add_subdirectory(dependencies/pbrt-parser EXCLUDE_FROM_ALL)
include_directories(dependencies/pbrt-parser/pbrtParser/include)

# OpenGL and X11 for the viewer
find_package(OpenGL REQUIRED)
if (LINUX)
        find_package(X11 REQUIRED)
        include_directories($(X11_INCLUDE_DIR))
	add_definitions(-DUSE_OPENGL)
elseif (WIN32)
	add_definitions(-DUSE_OPENGL)
endif()

# add the wurst sourcecode
file(GLOB_RECURSE wurst_src
    "src/*.h"
    "src/*.hpp"
    "src/*.cpp"
    "src/*.cc")
include_directories(src)
add_executable(wurstrenderer "${wurst_src}")

# link all external libraries
target_link_libraries(wurstrenderer embree assimp::assimp pbrtParser_static ${OPENGL_gl_LIBRARY})

# win32 vs other people
if (LINUX)
	target_link_libraries(wurstrenderer ${X11_LIBRARIES} GL stdc++fs)
endif ()

# for me
set(CMAKE_VERBOSE_MAKEFILE ON)
