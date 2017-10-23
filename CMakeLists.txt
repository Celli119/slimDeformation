cmake_minimum_required(VERSION 2.8.12)
project(slimDeformation)

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

FIND_PACKAGE(Matlab REQUIRED)
include_directories(${Matlab_INCLUDE_DIRS})

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)
option(LIBIGL_WITH_VIEWER ON)
option(LIBIGL_WITH_OPENGL ON)
add_subdirectory("${PROJECT_SOURCE_DIR}/libigl/shared/cmake" "libigl")


include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})

add_library(slim SHARED ${ADDITIONAL_SOURCES} ${SHARED_LIBRARIES} slimMatlab.cpp)
target_compile_definitions(slim PRIVATE ${LIBIGL_DEFINITIONS})
target_link_libraries(slim ${SHARED_LIBRARIES} ${MATLAB_LIBRARIES} ${LIBIGL_LIBRARIES} )


# rename the .so extension according to the OS
if (WIN32)
    set(mexExt ".mexw64")
elseif(APPLE)
    set(mexExt ".mexmaci64")
elseif(UNIX)
    set(mexExt ".mexa64")
endif()
set_target_properties(slim PROPERTIES SUFFIX ${mexExt})
