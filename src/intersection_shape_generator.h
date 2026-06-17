#pragma once
#include "types.h"
namespace isg {

/// 拓扑校验报告：包含错误与警告列表
struct ValidationReport {
    std::vector<std::string> errors, warnings;
    bool is_valid() const { return errors.empty(); }
};

/// 拓扑校验：检查车道组边界数量、连通关系转向类型与几何方位是否一致等
ValidationReport validateTopology(const IntersectionInput&);

/// 路口形态生成器：顶层入口，负责调度连通曲线生成、车道边线生成与路口面生成
class IntersectionShapeGenerator {
public:
    /// 配置参数
    struct Config {
        double sdf_cell_size = 0.2;       ///< SDF网格单元尺寸（米）
        double obstacle_buffer = 0.4;     ///< 障碍物缓冲距离（米）
        double kappa_max = 0.25;          ///< 最大允许曲率
        LBFGSConfig lbfgs;                ///< LBFGS优化器配置
        ConnectivityDirectionConfig connectivity_direction;  ///< 连通方向配置
    };

    IntersectionShapeGenerator();
    explicit IntersectionShapeGenerator(const Config& cfg);

    /// 生成路口形态：输入路口数据，输出连通曲线、车道边线与精细路口面
    /// @return true=成功 false=校验失败或异常
    bool generate(const IntersectionInput&, IntersectionOutput&);

    /// 获取最近一次校验报告
    const ValidationReport& lastReport() const { return report_; }

private:
    Config cfg_;
    ValidationReport report_;
};

}
