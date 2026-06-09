include(FetchContent)

# Catch2 v3 (header+source, MIT)
FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.4
    GIT_SHALLOW    TRUE)
FetchContent_MakeAvailable(Catch2)

# Helper to register a test executable
function(add_ig_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE
        intersection_gen_lib
        Catch2::Catch2WithMain
        #nlohmann_json::nlohmann_json
        )
    target_include_directories(${name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

enable_testing()

add_ig_test(test_bezier)
add_ig_test(test_sdf)
add_ig_test(test_optimizer)
add_ig_test(test_cluster)
add_ig_test(test_infeasibility)
add_ig_test(test_integration)
add_ig_test(test_scenarios)
