cmake_minimum_required(VERSION 3.1.0)
project(example)
set (CMAKE_CXX_STANDARD 20)

# add the Shape library
add_subdirectory(Shape)

# add the executable
add_executable(Example Example.cpp)

target_link_libraries(Example PUBLIC Shape)
