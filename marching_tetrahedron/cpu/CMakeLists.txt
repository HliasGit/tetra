file(GLOB_RECURSE TETRA_SOURCES src/*.c)
list(APPEND TETRA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../common/utils.c)

message(STATUS "TETRA_SOURCES: ${TETRA_SOURCES}")

add_library(TETRA ${TETRA_SOURCES})
target_include_directories(TETRA PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
)

if(VERBOSE)
    target_compile_definitions(TETRA PRIVATE VERBOSE)
endif()

target_link_libraries(TETRA PRIVATE COMMON)