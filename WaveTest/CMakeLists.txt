cmake_minimum_required(VERSION 3.0)
project(WaveTest)

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

set(WAVETEST_INCLUDE_FILES
    include/wave.h
   )
set(WAVETEST_SOURCE_FILES
    src/main.cpp
    src/wave.cpp
   )

set(SOURCE_FILES ${WAVETEST_INCLUDE_FILES} ${WAVETEST_SOURCE_FILES} ${PARAM_INC} ${PARAM_SRC})
add_executable(WaveTest ${SOURCE_FILES})
target_include_directories(WaveTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(WaveTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(WaveTest oops2d ${EXTRA_LIBS})
