find_package(OpenMP)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")


add_executable(intro_hello_world intro_hello_world.cpp)
target_link_libraries(intro_hello_world PRIVATE OpenMP::OpenMP_CXX)

add_executable(intro_sum intro_sum.cpp)
target_link_libraries(intro_sum PRIVATE OpenMP::OpenMP_CXX)

add_executable(intro_common_mem intro_common_mem.cpp)
target_link_libraries(intro_common_mem PRIVATE OpenMP::OpenMP_CXX)

add_executable(intro_common_mem_2 intro_common_mem_2.cpp)
target_link_libraries(intro_common_mem_2 PRIVATE OpenMP::OpenMP_CXX)

