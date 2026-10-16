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
    src/main.cpp
    src/simulation/simulation.cpp
    src/simulation/readParams.cpp
    src/structure/geometry.cpp
    src/structure/modedata.cpp
    src/structure/state.cpp
    src/structure/rungeStep.cpp
    src/fluid/channelSections.cpp
    src/fluid/flowModel.cpp
    src/fluid/separationSelection.cpp
    src/forces/forceCalculator.cpp
    src/forces/forceCalculatorFluid.cpp
    src/forces/forceCalculatorContact.cpp
    src/forces/modalProjector.cpp
    src/forces/surfaceLoad.cpp
)

set(SIMULATION_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/simulation
    ${CMAKE_CURRENT_SOURCE_DIR}/include/structure
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluid
    ${CMAKE_CURRENT_SOURCE_DIR}/include/forces
)

target_include_directories(simulation
    PRIVATE
        ${SIMULATION_INCLUDE_DIRS}
)

target_link_libraries(simulation
    PRIVATE
        OpenMP::OpenMP_CXX
)

if(PROFILE_CALCDIS)
    target_compile_definitions(simulation PRIVATE PROFILE_CALCDIS)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(test_separation_selection
        tests/test_separation_selection.cpp
        src/fluid/separationSelection.cpp
    )
    target_include_directories(test_separation_selection PRIVATE ${SIMULATION_INCLUDE_DIRS})
    add_test(NAME separation_selection COMMAND test_separation_selection)

    add_executable(test_mode_selection
        tests/test_mode_selection.cpp
        src/structure/modedata.cpp
    )
    target_include_directories(test_mode_selection PRIVATE ${SIMULATION_INCLUDE_DIRS})
    add_test(NAME mode_selection COMMAND test_mode_selection)

    add_executable(test_input_file_params
        tests/test_input_file_params.cpp
        src/simulation/readParams.cpp
    )
    target_include_directories(test_input_file_params PRIVATE ${SIMULATION_INCLUDE_DIRS})
    add_test(
        NAME input_file_params
        COMMAND test_input_file_params ${CMAKE_SOURCE_DIR}/input/param.txt
    )
endif()
