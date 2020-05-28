cmake_minimum_required(VERSION 3.0)
project(MeshTest)

if(NOT BUILD_TESTS)
  return()
endif()

# Generate the Parameters and ParamParser files.
set(PARAM_SRC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/meshparser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/meshparameters.cpp
   )
set(PARAM_INC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/meshparameters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/meshparser.h
   )
set(SETUP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mesh.json)

add_custom_command(
  OUTPUT ${PARAM_INC}
         ${PARAM_SRC}
  DEPENDS ${SETUP_SRC}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/genParams.py ${SETUP_SRC}
  COMMENT "Generating custom Parameters files"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM USES_TERMINAL
)

set(PARAMTEST_INCLUDE_FILES
   )
set(PARAMTEST_SOURCE_FILES
    src/main.cpp
   )

set(SOURCE_FILES ${PARAMTEST_INCLUDE_FILES} ${PARAMTEST_SOURCE_FILES} ${PARAM_INC} ${PARAM_SRC})
add_executable(MeshTest ${SOURCE_FILES})
target_include_directories(MeshTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(MeshTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(MeshTest oops2d ${EXTRA_LIBS})
