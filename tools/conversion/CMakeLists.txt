cmake_minimum_required(VERSION 3.10)
project(ConversionTool)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/../../include)

# Add the conversion tool executable.
add_executable(convert_serialized_graphene_data convert_serialized_graphene_data.cpp)

# Optionally, allow the user to set BLITZ_INCLUDE_DIR.
if(NOT DEFINED BLITZ_INCLUDE_DIR)
    set(BLITZ_INCLUDE_DIR "/usr/include" CACHE PATH "Path to Blitz headers")
endif()
include_directories(${BLITZ_INCLUDE_DIR})

# Check if Blitz is available by compiling an inline source snippet.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "
  #include <blitz/array.h>
  int main() {
    blitz::Array<int, 1> arr(10);
    return 0;
  }
  " BLITZ_WORKS
)

if(NOT BLITZ_WORKS)
  message(FATAL_ERROR "Could not compile test code including <blitz/array.h>. "
                      "Please check that Blitz is installed and BLITZ_INCLUDE_DIR is set correctly.")
else()
  message(STATUS "Blitz headers found and working.")
endif()

# Find Boost
find_package(Boost REQUIRED COMPONENTS serialization )
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(convert_serialized_graphene_data ${Boost_LIBRARIES})

# Add the compiler definition for boost serialization.
target_compile_definitions(convert_serialized_graphene_data  PRIVATE BZ_HAVE_BOOST_SERIALIZATION)
