cmake_minimum_required(VERSION 3.0)
project(medipix)

set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)


add_executable(medipix main.cpp Medipix.cpp)
target_link_libraries(medipix PRIVATE OpenMP::OpenMP_CXX)