cmake_minimum_required(VERSION 3.18.1)

if (USE_CUDA)
    project(pimc.e LANGUAGES CXX CUDA)
    add_definitions(-D USE_CUDA=${USE_CUDA})
else()
    project(pimc.e LANGUAGES CXX)
endif()

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

# Set c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_STANDARD_REQUIRED TRUE)
if (USE_CUDA)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 70)
    endif()
    message("CUDA architectures set to ${CMAKE_CUDA_ARCHITECTURES}.")
endif()

# Check for static build
if(STATIC)
    message("Static build specified, setting library suffixes to ${CMAKE_STATIC_LIBRARY_SUFFIX}.")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(BUILD_SHARED_LIBS OFF)
    list(APPEND CMAKE_EXE_LINKER_FLAGS "-static")
    list(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG "-static")
    set(Boost_USE_STATIC_LIBS   ON)
endif()

# Determine the dimension of space (default = 3)
if (NOT NDIM)
    set(NDIM 3)
endif()
add_definitions(-D NDIM=${NDIM})

# Check for GPU acceleration
if (GPU_BLOCK_SIZE)
    add_definitions(-D GPU_BLOCK_SIZE=${GPU_BLOCK_SIZE})
    if (NOT MAX_GPU_STREAMS)
        set(MAX_GPU_STREAMS 1)
    endif()
    add_definitions(-D MAX_GPU_STREAMS=${MAX_GPU_STREAMS})
endif()

# Set default build flags
if (NOT DEFAULT_CXX_FLAGS)
    set(DEFAULT_CXX_FLAGS "-Wall -fno-math-errno -O3")
    #set(DEFAULT_CXX_FLAGS "-march=native -Wall -fno-math-errno -O3")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEFAULT_CXX_FLAGS}")

#if (NOT EXTRA_CXX_FLAGS)
#    set(EXTRA_CXX_FLAGS "/Zc:__cplusplus")
#endif()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEFAULT_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")

# Determine executable name
set(exe pimc.e)
if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(exe pimcd.e)
elseif(CMAKE_BUILD_TYPE MATCHES PIGS)
    set(exe pigs.e)
elseif(CMAKE_BUILD_TYPE MATCHES PIGSDebug)
    set(exe pigsd.e)
endif()

# Find source files ( better to list explicitly https://stackoverflow.com/questions/1027247/specify-source-files-globally-with-glob )
file( GLOB PIMC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp )
if(USE_CUDA)
    file( GLOB PIMC_CUDA ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu )
    add_executable(${exe} ${PIMC_SRC} ${PIMC_CUDA})
    set_target_properties( ${exe} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    #if(APPLE)
    #  # We need to add the path to the driver (libcuda.dylib) as an rpath, 
    #  # so that the static cuda runtime can find it at runtime.
    #  set_property(TARGET ${exe} 
    #               PROPERTY
    #               BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    #endif()
else()
    add_executable(${exe} ${PIMC_SRC})
endif()

# Define headers for target
target_include_directories(${exe} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

# Set debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DPIMC_DEBUG")

# Add PIGS build mode
set(CMAKE_CXX_FLAGS_PIGS "${CMAKE_CXX_FLAGS} -DPIGS")
set(CMAKE_CXX_FLAGS_PIGS "${CMAKE_CXX_FLAGS_PIGS}" CACHE STRING
    "Flags used by the C++ compiler during PIGS builds."
    FORCE )
set(CMAKE_C_FLAGS_PIGS "" CACHE STRING
    "Flags used by the C compiler during PIGS builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_PIGS
    "${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING
    "Flags used for linking binaries during PIGS builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_PIGS
    "${CMAKE_SHARED_LINKER_FLAGS}" CACHE STRING
    "Flags used by the shared libraries linker during PIGS builds."
    FORCE )

# Add PIGSDEBUG build mode
set(CMAKE_CXX_FLAGS_PIGSDEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DPIGS")
set(CMAKE_CXX_FLAGS_PIGSDEBUG "${CMAKE_CXX_FLAGS_PIGSDEBUG}" CACHE STRING
    "Flags used by the C++ compiler during PIGS debug builds."
    FORCE )
set(CMAKE_C_FLAGS_PIGSDEBUG "${CMAKE_C_FLAGS_DEBUG}" CACHE STRING
    "Flags used by the C compiler during PIGS debug builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_PIGSDEBUG
    "${CMAKE_EXE_LINKER_FLAGS_DEBUG}" CACHE STRING
    "Flags used for linking binaries during PIGS debug builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_PIGSDEBUG
    "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" CACHE STRING
    "Flags used by the shared libraries linker during PIGS debug builds."
    FORCE )

mark_as_advanced(
    CMAKE_CXX_FLAGS_PIGS
    CMAKE_C_FLAGS_PIGS
    CMAKE_EXE_LINKER_FLAGS_PIGS
    CMAKE_SHARED_LINKER_FLAGS_PIGS
    CMAKE_CXX_FLAGS_PIGSDEBUG
    CMAKE_C_FLAGS_PIGSDEBUG
    CMAKE_EXE_LINKER_FLAGS_PIGSDEBUG
    CMAKE_SHARED_LINKER_FLAGS_PIGSDEBUG )

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE None
       CACHE STRING "Choose the type of build : None Debug Release PIGS PIGSDebug."
       FORCE)
endif()

if ((MSVC) AND (MSVC_VERSION GREATER_EQUAL 1914))
	target_compile_options(${exe} PUBLIC "/Zc:__cplusplus")
endif()

# Find Boost
find_package( Boost 1.71.0 REQUIRED COMPONENTS program_options serialization )

# Find Blitz
find_package( Blitz REQUIRED )
if (NOT Blitz_FOUND)
    message(FATAL_ERROR
        "Please follow blitz++ install instructions found at https://github.com/blitzpp/blitz.")
endif()

# Add include directories
include_directories( ${Boost_INCLUDE_DIRS} )
include_directories( ${Blitz_INCLUDE_DIRS} )

set(use_msvc_libs OFF)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if ("${CMAKE_CXX_SIMULATE_ID}" STREQUAL "MSVC")
        set(use_msvc_libs ON)
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(use_msvc_libs ON)
endif ()

# Link libraries
target_link_libraries (${exe} ${Boost_LIBRARIES} $<$<BOOL:${use_msvc_libs}>:bcrypt>)
if((CMAKE_BUILD_TYPE MATCHES Debug) OR (CMAKE_BUILD_TYPE MATCHES PIGSDebug))
    target_link_libraries (${exe} ${Blitz_LIBRARIES})
endif()

# Link filesystem library -lstdc++fs for old compilers
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0))
    message("Linking filesystem libraries -lstdc++fs for older compilers")
    target_link_libraries(${exe} stdc++fs)
endif()

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${exe} EXPORT ${exe}Config
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.
install(EXPORT ${exe}Config DESTINATION share/${exe}/cmake)

# This makes the project importable from the build directory
export(TARGETS ${exe} FILE ${exe}Config.cmake)

# Unit tests
# We need to make some of these
enable_testing()
