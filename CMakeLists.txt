cmake_minimum_required(VERSION 3.0)
project(medipix)

set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)

include_directories(PkgConfig::FFTW)
include_directories(include)

add_library(medipix src/Medipix.cpp src/helper.cpp src/MedipixSPM.cpp src/MedipixCSM.cpp )
target_link_libraries(medipix PRIVATE OpenMP::OpenMP_CXX PkgConfig::FFTW)

add_subdirectory(tests)

add_subdirectory(example)