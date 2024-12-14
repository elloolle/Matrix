cmake_minimum_required(VERSION 3.28)
project(Matrix)

set(CMAKE_CXX_STANDARD 23)

add_executable(Matrix main.cpp MatrixClass.h
        Residue.h Rational.h matrix.h)
add_compile_options(-std=c++23)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=undefined")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=undefined")
