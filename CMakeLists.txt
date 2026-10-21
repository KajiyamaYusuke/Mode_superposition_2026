cmake_minimum_required(VERSION 3.16)

# project() より前に置く方が安全
set(CMAKE_CXX_COMPILER icpx)

project(Simulation LANGUAGES CXX)

option(PROFILE_CALCDIS "Enable calcDis profiling output" OFF)
option(ENABLE_KANAYA_DIAGNOSTICS "Build opt-in read-only force/work diagnostics" OFF)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/include/diagnostics
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

if(ENABLE_KANAYA_DIAGNOSTICS)
    target_sources(simulation PRIVATE src/diagnostics/KanayaDiagnostics.cpp)
    target_compile_definitions(simulation PRIVATE ENABLE_KANAYA_DIAGNOSTICS)
endif()

include(CTest)
if(BUILD_TESTING)
    if(ENABLE_KANAYA_DIAGNOSTICS)
        get_target_property(KANAYA_TEST_SOURCES simulation SOURCES)
        list(REMOVE_ITEM KANAYA_TEST_SOURCES src/main.cpp)
        add_executable(test_kanaya_diagnostics tests/test_kanaya_diagnostics.cpp ${KANAYA_TEST_SOURCES})
        target_include_directories(test_kanaya_diagnostics PRIVATE ${SIMULATION_INCLUDE_DIRS})
        target_compile_definitions(test_kanaya_diagnostics PRIVATE ENABLE_KANAYA_DIAGNOSTICS)
        target_link_libraries(test_kanaya_diagnostics PRIVATE OpenMP::OpenMP_CXX)
        add_test(NAME kanaya_diagnostics COMMAND test_kanaya_diagnostics
            ${CMAKE_SOURCE_DIR}/output/kanaya_diagnostic/validation/unit)
        set_tests_properties(kanaya_diagnostics PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")
    else()
        get_target_property(KANAYA_TEST_SOURCES simulation SOURCES)
        list(REMOVE_ITEM KANAYA_TEST_SOURCES src/main.cpp)
        add_executable(kanaya_smoke_driver tests/test_kanaya_diagnostics.cpp ${KANAYA_TEST_SOURCES})
        target_include_directories(kanaya_smoke_driver PRIVATE ${SIMULATION_INCLUDE_DIRS})
        target_link_libraries(kanaya_smoke_driver PRIVATE OpenMP::OpenMP_CXX)
    endif()
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
