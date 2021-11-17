cmake_minimum_required(VERSION 3.10)

# set project name and version
project(Multigrid VERSION 0.1.0)

# specify the language
enable_language(CXX)

# specify compliers for C and C++ 
set(CMAKE_C_COMPILER=gcc)
set(CMAKE_CXX_COMPILER=g++)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(Multigrid Prolong_restrict.cpp)
