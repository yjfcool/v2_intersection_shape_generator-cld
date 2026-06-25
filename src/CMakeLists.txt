add_library(intersection_gen_lib STATIC
        intersection_shape_generator.cpp

        utils/clipper.cpp
        utils/quadtree.cpp

        curve/bezier.cpp
        curve/hermite_init.cpp
        curve/curve_utils.cpp

        optimizer/sdf_field.cpp
        optimizer/penalty_cost.cpp
        optimizer/lbfgs_solver.cpp

        constraints/infeasibility_detector.cpp
        constraints/cluster_order.cpp

        generator/connectivity_generator.cpp
        generator/polygon_builder.h
        generator/edge_line_generator.h
)

target_include_directories(intersection_gen_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

if (DEPENDS_EXISTS STREQUAL "TRUE")
    target_link_libraries(intersection_gen_lib PUBLIC)
else ()
    target_link_libraries(intersection_gen_lib PUBLIC Eigen3::Eigen)
endif ()

target_compile_options(intersection_gen_lib PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra>
        $<$<CXX_COMPILER_ID:MSVC>:/W3>
)

# 性能优化: 不在CMakeLists中硬编码-O2, 由CMAKE_BUILD_TYPE控制
# - Debug: -g -Og (保留Debug快速分支,不定义NDEBUG)
# - Release: -O3 -DNDEBUG
# - RelWithDebInfo: -O2 -g -DNDEBUG
target_compile_options(intersection_gen_lib PRIVATE
        $<$<CONFIG:Debug>:-g -Og -DNDEBUG>
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
        $<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
)

# 添加测试用例应用程序
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
    add_executable(intersection_gen main.cpp)
    target_link_libraries(intersection_gen PRIVATE intersection_gen_lib)
    # 设置测试应用程序的编译选项
    target_compile_options(intersection_gen PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra>
            $<$<CXX_COMPILER_ID:MSVC>:/W3>
    )
endif ()
