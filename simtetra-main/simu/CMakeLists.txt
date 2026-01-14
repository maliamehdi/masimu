cmake_minimum_required(VERSION 3.16)
project(simTetra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED ui_all vis_all gdml)

if(TARGET Geant4::geant4)
  set(GEANT4_LIBS Geant4::geant4)
else()
  include(${Geant4_USE_FILE})
  set(GEANT4_LIBS ${Geant4_LIBRARIES})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIRS})

file(GLOB SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc
)

set(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/simTetra.cc
  ${SRC_FILES}
)

add_executable(simTetra ${SOURCES})
target_link_libraries(simTetra PRIVATE ${GEANT4_LIBS})

install(TARGETS simTetra RUNTIME DESTINATION bin)
