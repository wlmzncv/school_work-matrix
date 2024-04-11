cmake_minimum_required(VERSION 3.24.1)
# 生成Makefile文件到build目录:
# cmake -G "MinGW Makefiles" -B build .
# 利用Makefile文件生成exe:
# cmake --build ./build
project(Matrix)

aux_source_directory(. DIR_SRCS)

add_executable(Matrix ${DIR_SRCS})