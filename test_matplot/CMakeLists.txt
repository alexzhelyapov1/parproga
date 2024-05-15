cmake_minimum_required(VERSION 3.16)
project("PARPROGA 1")
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(matplotplusplus)

add_executable(main main.cpp)
target_link_libraries(main matplot)