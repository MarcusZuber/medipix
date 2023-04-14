cmake_minimum_required(VERSION 3.0)
project(medipix)

set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

include_directories(include)

add_library(medipix src/Medipix.cpp src/helper.cpp src/MedipixSPM.cpp src/MedipixCSM.cpp )
target_link_libraries(medipix PRIVATE OpenMP::OpenMP_CXX)

add_subdirectory(tests)

add_subdirectory(example)