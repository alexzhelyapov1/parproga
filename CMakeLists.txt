cmake_minimum_required(VERSION 3.16)
project("PARPROGA 1")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

SET(MYCOMPILE_FLAGS "-pthread")
SET(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS} ${MYCOMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MYCOMPILE_FLAGS}")

find_package(MPI)

file(GLOB SRC "code/src/*")

add_subdirectory(7_sem)

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


add_custom_target(1_lab ALL
  COMMAND mpic++ ${CMAKE_CURRENT_SOURCE_DIR}/src/1_lab.cpp -o ${CMAKE_CURRENT_BINARY_DIR}/1_lab
) # mpiexec -n 4 ./1_lab

add_executable(1_lab_thread src/1_lab_1_thread.cpp)

# add_executable(1_lab_new src/1_lab_new.cpp)
add_custom_target(1_lab_new ALL
  COMMAND mpic++ ${CMAKE_CURRENT_SOURCE_DIR}/src/1_lab_new.cpp -o ${CMAKE_CURRENT_BINARY_DIR}/1_lab_new
) # mpiexec -n 2 ./1_lab_new


# add_custom_command(TARGET MAIN_PROG POST_BUILD
#   COMMAND mpiexec -n 2 ${CMAKE_CURRENT_BINARY_DIR}/main
#   COMMENT "Запуск MAIN_PROG..."
# )

# add_executable(test src/test.cpp)
add_custom_target(mpi_test ALL
  COMMAND mpic++ ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi_test.cpp -o ${CMAKE_CURRENT_BINARY_DIR}/mpi_test
) # mpiexec -n 4 ./mpi_test


add_custom_target(
    compare_solutions
    COMMAND diff -q /home/alex/parproga/gen/solution_new.csv /home/alex/parproga/gen/solution_new_new.csv
    VERBATIM
)

# add_executable(test_matplot src/test_matplot.cpp)
# target_link_libraries(test_matplot PUBLIC matplot)
