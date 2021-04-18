cmake_minimum_required(VERSION 3.20.1)

project(pimc.e)

# Set c++ standard
set(CMAKE_CXX_STANDARD 17)

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
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -fno-math-errno -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fno-math-errno -O3")

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
file(GLOB PIMC_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    )
add_executable(${exe} ${PIMC_SRC})

# Define headers for target
target_include_directories(${exe} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

# Set debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D PIMC_DEBUG")

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

# Find Boost
find_package( Boost 1.74.0 REQUIRED COMPONENTS program_options filesystem serialization )

# Find Blitz
find_package( Blitz REQUIRED )
if (NOT Blitz_FOUND)
    message(FATAL_ERROR
        "Please follow blitz++ install instructions found at https://github.com/blitzpp/blitz.")
endif()

# Add include directories
include_directories( ${Boost_INCLUDE_DIRS} )
include_directories( ${Blitz_INCLUDE_DIRS} )

# Link libraries
target_link_libraries (${exe} ${Boost_LIBRARIES} )
if((CMAKE_BUILD_TYPE MATCHES Debug) OR (CMAKE_BUILD_TYPE MATCHES PIGSDebug))
    target_link_libraries (${exe} ${Blitz_LIBRARIES} )
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
