cmake_minimum_required(VERSION 3.0)
project(PolarTest)

if(NOT BUILD_TESTS)
  return()
endif()

# Generate the Parameters and ParamParser files.
set(PARAM_SRC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/waveparser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/waveparameters.cpp
   )
set(PARAM_INC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/waveparameters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/waveparser.h
   )
set(SETUP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/scripts/wave.json)

add_custom_command(
  OUTPUT ${PARAM_INC}
         ${PARAM_SRC}
  DEPENDS ${SETUP_SRC}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/genParams.py ${SETUP_SRC}
  COMMENT "Generating custom Parameters files"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM USES_TERMINAL
)

set(POLARTEST_INCLUDE_FILES
    include/wave.h
   )
set(POLARTEST_SOURCE_FILES
    src/main.cpp
    src/wave.cpp
   )

set(SOURCE_FILES ${POLARTEST_INCLUDE_FILES} ${POLARTEST_SOURCE_FILES} ${PARAM_INC} ${PARAM_SRC})
add_executable(PolarTest ${SOURCE_FILES})
target_include_directories(PolarTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(PolarTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(PolarTest oops2d ${EXTRA_LIBS})
