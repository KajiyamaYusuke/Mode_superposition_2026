cmake_minimum_required(VERSION 3.16)
project(Simulation LANGUAGES CXX)

set(CMAKE_CXX_COMPILER icpx)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ipo -xHost")
#set(CMAKE_BUILD_TYPE Debug)

# include ディレクトリを追加
include_directories(include)

add_executable(simulation
    src/forceCalculator.cpp
    src/geometry.cpp
    src/modedata.cpp
    src/readParams.cpp
    src/rungeStep.cpp
    src/simulation.cpp
    src/state.cpp
    src/main.cpp          # main を含むファイル
)



