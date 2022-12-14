cmake_minimum_required(VERSION 3.13.0)

project(superlu C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../lib)

file(GLOB SUPERLU_SRC "*.c")

add_compile_definitions(NDEBUG)
add_compile_definitions(PRNTlevel=0)
add_compile_definitions(Add_)

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
    link_libraries(libopenblas)
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    link_libraries(pthread gfortran quadmath openblas)
endif ()

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SUPERLU_SRC})

message(STATUS "SUPERLU C_FLAGS: ${CMAKE_C_FLAGS}")
