cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(RoutingFramework CXX)

# Directories containing .cc source files.
set(SOURCE_DIRECTORIES Launchers RawData)

# Determine the language standard.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Flags when building for the Devel configuration.
if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
  set(CMAKE_CXX_FLAGS_DEVEL -O3)
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Devel)
endif()

# Enable the compiler to use extended instructions in generated code.
if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
option(DISABLE_AVX "Disable use of instructions in the AVX extended instruction set.")
if(DISABLE_AVX)
  if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-avx")
  endif()
endif()

# Request warnings.
if(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wfatal-errors -Wpedantic -pedantic-errors\
      -Wall -Wextra -ftemplate-backtrace-limit=1 -Wno-unknown-pragmas -Wno-sign-compare")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic-errors -ferror-limit=1\
      -ftemplate-backtrace-limit=1 -Wno-sign-compare")
endif()

include_directories(${CMAKE_SOURCE_DIR})

find_package(OpenMP)
find_package(CGAL REQUIRED)
find_library(cairo_LIBRARY cairo)
find_library(PROJ_LIBRARY proj)
find_library(RoutingKit_LIBRARY routingkit)

# Clang does not support the option -frounding-math, thus remove it from the CGAL compile options.
if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set_target_properties(CGAL::CGAL PROPERTIES INTERFACE_COMPILE_OPTIONS "")
endif()

foreach(SRC_DIR ${SOURCE_DIRECTORIES})
  add_subdirectory(${SRC_DIR})
endforeach()
