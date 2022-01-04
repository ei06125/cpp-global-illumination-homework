cmake_minimum_required(VERSION 3.21)
project(smallpt)

find_package(OpenMP REQUIRED)

add_executable(smallpt)
target_sources(smallpt PRIVATE src/smallpt.cpp)
# target_link_libraries(smallpt PRIVATE OpenMP)
target_compile_definitions(smallpt PRIVATE -fopenmp)
