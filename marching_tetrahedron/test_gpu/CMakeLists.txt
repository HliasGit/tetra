# Macro to add tests
macro(add_test_executable_c name)
    add_executable(${name} ${name}.cu)
    target_link_libraries(${name} PUBLIC TETRA_GPU m)
    target_include_directories(${name} PUBLIC 
        ${CMAKE_SOURCE_DIR}/gpu/inc
        ${CMAKE_SOURCE_DIR}/common
    )

    set_target_properties(${name} PROPERTIES CUDA_ARCHITECTURES "86")

    if(VERBOSE)
        target_compile_definitions(${name} PUBLIC VERBOSE)
    endif()
    message(STATUS "Including test ${name}")
endmacro()

# Test executables
add_test_executable_c(test)
add_test_executable_c(test_wo_preprocessing)
add_test_executable_c(test_optimized)
add_test_executable_c(test_tetra_based)
add_test_executable_c(test_flags_noatom)
add_test_executable_c(test_mixed)