cmake_minimum_required(VERSION 3.0)

# Set the project name
project(OOPS2D VERSION 1.0)

# Specify the C++ Standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Require Python
find_package(PythonInterp 3.0 REQUIRED)
find_package(MPI REQUIRED)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

# Set a default build type.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to Debug as none was specified.")
  set(CMAKE_BUILD_TYPE "Debug" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(BUILD_TESTS "Compile OOPS unit tests" OFF)

#if(NOT USE_AVX)
#  set(USE_AVX "None" CACHE
#      STRING "Compile with AVX. This usually speeds up code, but it may not be available on all machines." FORCE)
#  set_property(CACHE USE_AVX PROPERTY STRINGS
#    "None" "AVX" "AVX2" "AVX-512")
#endif()

option(USE_INTRINSICS "On GNU-like and Intel compilers, use processor-specific instruction sets to improve performance." OFF)

# Configure a header file that will pass on our settings.
configure_file(${PROJECT_SOURCE_DIR}/include/oopsconfig.h.in ${PROJECT_SOURCE_DIR}/include/oopsconfig.h)

# A macro to do a recursive search for subdirectories with CMakeLists.txt files.
macro(SUBDIRLIST result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child} AND EXISTS ${curdir}/${child}/CMakeLists.txt)
      list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

set(EXTRA_LIBS
    ${MPI_LIBRARIES}
   )

if(USE_INTRINSICS)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-march=native)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if(UNIX)
      add_compile_options(-xHost)
    else()
      add_compile_options(/QxHost)
    endif()
  else()
    message(STATUS "USE_INTRINSICS is not set up for your compiler and therefore does nothing.")
  endif()
endif()

# Add the AVX flag for the different compilers.
#if(NOT "${USE_AVX}" STREQUAL "None")
#  if("${USE_AVX}" STREQUAL "AVX")
#    if(MSVC)
#      add_compile_options(/arch:AVX)
#    else()
#      add_compile_options(-mavx)
#    endif()
#  endif()
#  if("${USE_AVX}" STREQUAL "AVX2")
#    if(MSVC)
#      add_compile_options(/arch:AVX2)
#    else()
#      add_compile_options(-mavx2)
#    endif()
#  endif()
#  if("${USE_AVX}" STREQUAL "AVX512")
#    if(MSVC)
#      add_compile_options(/arch:AVX512)
#    else()
#      add_compile_options(-mavx512f)
#    endif()
#  endif()
#endif()

set(OOPS_INCLUDE_FILES
    include/mpicommunicator.h
    include/types.h
    include/parameters.h
    include/paramreader.h
    include/paramparser.h
    include/grid.h
    include/domain.h
    include/ode.h
    include/odedata.h
    include/solverdata.h
    include/fieldmap.h
    include/solver.h
    include/output.h
    include/rk4.h
    include/rkck.h
    include/base64.h
   )

set(OOPS_SOURCE_FILES
    src/mpicommunicator.cpp
    src/paramreader.cpp
    src/grid.cpp
    src/domain.cpp
    src/fieldmap.cpp
    src/ode.cpp
    src/odedata.cpp
    src/solverdata.cpp
    src/output.cpp
    src/rk4.cpp
    src/rkck.cpp
    src/base64.cpp
   )

add_library(oops2d ${OOPS_INCLUDE_FILES} ${OOPS_SOURCE_FILES})
target_include_directories(oops2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(oops2d PUBLIC ${MPI_INCLUDE_PATH})

# Automatically add subdirectories containing CMakeLists.txt files and compile them.
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
foreach(subdir ${SUBDIRS})
  add_subdirectory(${subdir})
endforeach()
