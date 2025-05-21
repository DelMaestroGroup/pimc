cmake_minimum_required(VERSION 3.18.1)

project(pimc.e LANGUAGES CXX)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

# Set c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_STANDARD_REQUIRED TRUE)

include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-std=c++17")

# --------- Filesystem Header Test ---------
check_cxx_source_compiles(
  "
  #include <filesystem>
  int main() {
    std::filesystem::path p;
    return 0;
  }
  " HAVE_STD_FILESYSTEM)

if (HAVE_STD_FILESYSTEM)
  message(STATUS "Found <filesystem>")
  add_compile_definitions(USE_STD_FILESYSTEM)
else()
  check_cxx_source_compiles(
    "
    #include <experimental/filesystem>
    int main() {
      std::experimental::filesystem::path p;
      return 0;
    }
    " HAVE_EXPERIMENTAL_FILESYSTEM)
  if (HAVE_EXPERIMENTAL_FILESYSTEM)
    message(STATUS "Found <experimental/filesystem>")
    add_compile_definitions(USE_EXPERIMENTAL_FILESYSTEM)
  else()
    message(FATAL_ERROR "Neither <filesystem> nor <experimental/filesystem> are available")
  endif()
endif()
# --------- End Filesystem Header Test ---------

# --------- MDSpan Header Test ---------
# === 1) Check for standard <mdspan> ===
check_cxx_source_compiles(
  "
  #include <mdspan>
  int main() {
    std::mdspan<float, std::dextents<size_t,2>> m(nullptr, 3, 4);
    (void)m;
    return 0;
  }
  "
  HAVE_STD_MDSPAN
)

if(HAVE_STD_MDSPAN)
  message(STATUS "Found <mdspan>")
  add_compile_definitions(USE_STD_MDSPAN)
else()
  # === 2) Check for <experimental/mdspan> ===
  check_cxx_source_compiles(
    "
    #include <experimental/mdspan>
    int main() {
      std::experimental::mdspan<float, std::experimental::dextents<size_t,2>> m(nullptr, 3, 4);
      (void)m;
      return 0;
    }
    "
    HAVE_EXPERIMENTAL_MDSPAN
  )

  if(HAVE_EXPERIMENTAL_MDSPAN)
    message(STATUS "Found <experimental/mdspan>")
    add_compile_definitions(USE_EXPERIMENTAL_MDSPAN)
  else()
    # === 3) Download the single-header reference implementation ===

    # (A) Define where to place the reference header
    set(REFERENCE_MDSPAN_NAME "mdspan.h")
    set(REFERENCE_MDSPAN_BINARY "${CMAKE_BINARY_DIR}/mdspan_downloaded.hpp")
    set(REFERENCE_MDSPAN_INSTALL "${CMAKE_CURRENT_SOURCE_DIR}/include/${REFERENCE_MDSPAN_NAME}")

    # (B) Download using CMake's file(DOWNLOAD) command.
    file(DOWNLOAD
      "https://raw.githubusercontent.com/kokkos/mdspan/single-header/mdspan.hpp"
      "${REFERENCE_MDSPAN_BINARY}"
      EXPECTED_MD5 "72604e49ef586d353c51257eb30ff027"  # optional hash check
      STATUS DOWNLOAD_STATUS
    )

    list(GET DOWNLOAD_STATUS 0 DOWNLOAD_STATUS_CODE)
    list(GET DOWNLOAD_STATUS 1 DOWNLOAD_STATUS_STRING)

    if(NOT ${DOWNLOAD_STATUS_CODE} EQUAL 0)
      message(FATAL_ERROR "mdspan.hpp download failed: ${DOWNLOAD_STATUS_STRING}")
    endif()

    # (C) Create the local include directory if it doesn't exist
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include")

    # (D) Rename (copy) the downloaded file to mdspan.h in our include dir
    file(COPY
      "${REFERENCE_MDSPAN_BINARY}"
      DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
    )

    file(RENAME
      "${CMAKE_CURRENT_SOURCE_DIR}/include/mdspan_downloaded.hpp"
      "${REFERENCE_MDSPAN_INSTALL}"
    )

    # (E) Test if the newly downloaded file compiles
    check_cxx_source_compiles(
      "
      #include \"${REFERENCE_MDSPAN_INSTALL}\"
      int main() {
        std::mdspan<float, std::dextents<size_t,2>> m(nullptr, 3, 4);
        (void)m;
        return 0;
      }
      "
      HAVE_REFERENCE_MDSPAN
    )

    if(HAVE_REFERENCE_MDSPAN)
      message(STATUS "Downloaded and verified the mdspan reference implementation.")
      add_compile_definitions(USE_REFERENCE_MDSPAN)
    else()
      message(FATAL_ERROR
        "Neither <mdspan>, <experimental/mdspan>, nor the downloaded reference "
        "implementation compiled successfully.")
    endif()
  endif()
endif()

# --------- End MDSpan Header Test ---------

# Enable GPU support
set(GPU_COMPILE_FLAGS "")
set(GPU_BACKEND "none" CACHE STRING "Enable gpu accleration: cuda, hip, sycl, none (default: none)")
if (NOT ${GPU_BACKEND} STREQUAL "none")
    add_definitions(-D USE_GPU=1)
    if (${GPU_BACKEND} STREQUAL "cuda")
        add_definitions(-D USE_CUDA=1)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES 80)
        endif()
        message("CUDA architectures set to ${CMAKE_CUDA_ARCHITECTURES}.")
        enable_language(CUDA)
        include_directories("${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    elseif(${GPU_BACKEND} STREQUAL "hip")
        add_definitions(-D USE_HIP=1)
    elseif(${GPU_BACKEND} STREQUAL "sycl")
        add_definitions(-D USE_SYCL=1)
        if (NOT SYCL_FLAGS)
            list(APPEND GPU_COMPILE_FLAGS
                "-fsycl"
            )
        else()
            list(APPEND GPU_COMPILE_FLAGS
                ${SYCL_FLAGS}
            )
        endif()
    else()
        message(FATAL_ERROR "Unknown GPU implementation '${GPU_BACKEND}',"
            " please select from 'cuda, hip, sycl, none'.")
    endif()

    # Set number of threads per block
    set(GPU_BLOCK_SIZE "1024" CACHE STRING "Set gpu block size (default: 1024)")
    add_definitions(-D GPU_BLOCK_SIZE=${GPU_BLOCK_SIZE})
    
    # Set sub-group size (warpsize/wavefront/SIMD lanes/etc.)
    if(${GPU_BACKEND} STREQUAL "hip")
        set(SUB_GROUP_SIZE "64" CACHE STRING "Set sub-group size (default: 32)")
        message("HIP detected. Set SUB_GROUP_SIZE to ${SUB_GROUP_SIZE}.")
    else()
        set(SUB_GROUP_SIZE "32" CACHE STRING "Set sub-group size (default: 32)")
    endif()
    add_definitions(-D SUB_GROUP_SIZE=${SUB_GROUP_SIZE})
    
    # Set number of GPU streams
    set(MAX_GPU_STREAMS "1" CACHE STRING "Set number of gpu streams (default: 1)")
    add_definitions(-D MAX_GPU_STREAMS=${MAX_GPU_STREAMS})
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

# Determine if we want to compile for boltzmannons
if (BOLTZMANNONS)
    add_compile_definitions(BOLTZMANNONS)
endif()

# Set default build flags
if (NOT DEFAULT_CXX_FLAGS)
    set(DEFAULT_CXX_FLAGS "-Wall -fno-math-errno -O3")
    #set(DEFAULT_CXX_FLAGS "-march=native -Wall -fno-math-errno -O3")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEFAULT_CXX_FLAGS}")

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
if (${GPU_BACKEND} STREQUAL "cuda")
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

# Find Boost
find_package( Boost 1.71.0 REQUIRED COMPONENTS program_options serialization )

# Add include directories
include_directories( ${Boost_INCLUDE_DIRS} )

# Add compile flags
target_compile_options(${exe} PRIVATE ${GPU_COMPILE_FLAGS})

# Link libraries
target_link_libraries (${exe} ${Boost_LIBRARIES} )

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
