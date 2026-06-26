cmake_minimum_required(VERSION 3.16)

# project() より前に置く方が安全
set(CMAKE_CXX_COMPILER icpx)

project(Simulation LANGUAGES CXX)

option(PROFILE_CALCDIS "Enable calcDis profiling output" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# single-config generator 用
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Release 最適化
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)

add_executable(simulation
    src/forceCalculator.cpp
    src/geometry.cpp
    src/modedata.cpp
    src/readParams.cpp
    src/rungeStep.cpp
    src/simulation.cpp
    src/state.cpp
    src/main.cpp
)

target_include_directories(simulation
    PRIVATE
        include
)

target_link_libraries(simulation
    PRIVATE
        OpenMP::OpenMP_CXX
)

if(PROFILE_CALCDIS)
    target_compile_definitions(simulation PRIVATE PROFILE_CALCDIS)
endif()