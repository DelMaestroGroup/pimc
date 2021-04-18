cmake_minimum_required(VERSION 3.5.1)
project(pimc.e LANGUAGES CXX)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

# Targets to develop
enable_testing()
add_subdirectory(pimc)
