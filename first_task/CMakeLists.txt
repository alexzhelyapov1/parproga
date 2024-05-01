cmake_minimum_required(VERSION 3.16)
project("PARPROGA 1")
set(CMAKE_CXX_STANDARD 17)
SET(MYCOMPILE_FLAGS "-pthread")

SET(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS} ${MYCOMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MYCOMPILE_FLAGS}")

find_package(MPI)

file(GLOB SRC "code/src/*")



add_custom_target(PI ALL
  COMMAND mpic++ ${CMAKE_CURRENT_SOURCE_DIR}/src/pi.cpp -o ${CMAKE_CURRENT_BINARY_DIR}/pi
) # mpiexec -n 16 ./pi


add_custom_target(COMMUNICATION_TIME ALL
  COMMAND mpic++ ${CMAKE_CURRENT_SOURCE_DIR}/src/communication_time.cpp -o ${CMAKE_CURRENT_BINARY_DIR}/com_time
) # mpiexec -n 16 ./com_time


# add_custom_target(INTEGRAL ALL
#   COMMAND mpic++ ${CMAKE_CURRENT_SOURCE_DIR}/src/integral.cpp -o ${CMAKE_CURRENT_BINARY_DIR}/integral
# ) # mpiexec -n 16 ./integral

add_executable(integral src/integral.cpp)
# target_compile_options()


# add_custom_command(TARGET MAIN_PROG POST_BUILD
#   COMMAND mpiexec -n 2 ${CMAKE_CURRENT_BINARY_DIR}/main
#   COMMENT "Запуск MAIN_PROG..."
# )

add_executable(test src/trash.cpp)