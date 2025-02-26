cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_C_COMPILER clang-omp)
# set(CMAKE_CPP_COMPILER clang-omp++)

# message(STATUS "set env var")
# set(ENV{CC} "/usr/local/opt/llvm/bin/clang")
# set(ENV{CXX} "/usr/local/opt/llvm/bin/clang++")
# set(ENV{LDFLAGS} "-L/usr/local/opt/llvm/lib")
# set(ENV{CPPFLAGS} "-I/usr/local/opt/llvm/include")

project(bfield_modelling VERSION 1.0)

file(GLOB_RECURSE SRC_FILES src/*.cpp)
add_executable(bfield_modelling ${SRC_FILES})

target_include_directories(bfield_modelling PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(OpenMP REQUIRED)
target_link_libraries(bfield_modelling PUBLIC OpenMP::OpenMP_CXX)

# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()