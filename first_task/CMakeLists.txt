cmake_minimum_required(VERSION 3.16)
project("PARPROGA 1")
set(CMAKE_CXX_STANDARD 23)
SET(MYCOMPILE_FLAGS "-pthread")

SET(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS} ${MYCOMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MYCOMPILE_FLAGS}")

find_package(MPI)

file(GLOB SRC "code/src/*")

add_executable(parproga_1 ${SRC})

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

target_link_libraries(parproga_1 ${MPI_C_LIBRARIES})
target_include_directories(parproga_1 PUBLIC "code/include")