#pragma once
#include "types.h"
#include "curve/bezier.h"
#include <vector>
#include <memory>

namespace isg {

struct BoundingBox2d;

/// C++11兼容的make_unique
template<typename T, typename... Args>
std::unique_ptr<T> make_unique_cpp11(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/// 四叉树空间索引：用于加速Bezier曲线的范围查询
class QuadTree {
public:
    /// 四叉树节点
    struct Node {
        BoundingBox2d bounds;
        std::vector<std::pair<BezierCurve, std::string>> curves;  ///< (曲线, ID)
        std::unique_ptr<Node> children[4];  ///< NW, NE, SW, SE 四个子节点
        bool divided = false;

        Node(const BoundingBox2d& b) : bounds(b) {
            children[0] = nullptr;
            children[1] = nullptr;
            children[2] = nullptr;
            children[3] = nullptr;
        }
    };

private:
    std::unique_ptr<Node> root_;
    size_t capacity_;     ///< 节点容量，超过则细分
    size_t depth_limit_;  ///< 最大深度

public:
    QuadTree(const BoundingBox2d& bounds, size_t capacity = 8, size_t depth_limit = 8)
        : capacity_(capacity), depth_limit_(depth_limit) {
        root_ = make_unique_cpp11<Node>(bounds);
    }

    /// 插入曲线
    void insert(const BezierCurve& curve, const std::string& id);

    /// 范围查询：返回与查询包围盒相交的所有曲线
    std::vector<std::pair<BezierCurve, std::string>> queryRange(const BoundingBox2d& range) const;

    /// 清空四叉树
    void clear();

private:
    void insertRecursive(Node* node, const BezierCurve& curve, const std::string& id, size_t current_depth);
    void queryRecursive(const Node* node, const BoundingBox2d& range,
                       std::vector<std::pair<BezierCurve, std::string>>& results) const;
    void subdivide(Node* node);
    bool intersects(const BoundingBox2d& box, const BezierCurve& curve) const;
    BoundingBox2d getCurveBounds(const BezierCurve& curve) const;
};

}
