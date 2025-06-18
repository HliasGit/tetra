# Macro to add tests
macro(add_test_executable_c name)
    add_executable(${name} ${name}.c)
    target_link_libraries(${name} PUBLIC TETRA m)
    target_include_directories(${name} PUBLIC 
        ${CMAKE_SOURCE_DIR}/cpu/inc
        ${CMAKE_SOURCE_DIR}/common
    )

    if(VERBOSE)
        target_compile_definitions(${name} INTERFACE VERBOSE)
    endif()
    message(STATUS "Including test ${name}")
endmacro()

# Macro to add tests
macro(add_test_executable_cpp name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PUBLIC TETRA m)
    target_include_directories(${name} PUBLIC 
        ${CMAKE_SOURCE_DIR}/cpu/inc
    )

    if(VERBOSE)
        target_compile_definitions(${name} INTERFACE VERBOSE)
    endif()
    message(STATUS "Including test ${name}")
endmacro()

# Test executables
add_test_executable_c(main)
add_test_executable_c(non_unique)
add_test_executable_c(generator)

add_test_executable_cpp(normals)
add_test_executable_cpp(euler_characteristic)