#include "quadtree.h"
#include "utils.h"
#include <algorithm>
#include <cmath>

void QuadTree::subdivide(Node* node) {
    double half_width = node->bounds.width() / 2.0;
    double half_height = node->bounds.height() / 2.0;
    Vec2d center = (node->bounds.min_pt + node->bounds.max_pt) / 2.0;

    // Create four quadrants: NW, NE, SW, SE
    BoundingBox2d nw_bounds;
    nw_bounds.min_pt = Vec2d(node->bounds.min_pt.x(), center.y());
    nw_bounds.max_pt = Vec2d(center.x(), node->bounds.max_pt.y());

    BoundingBox2d ne_bounds;
    ne_bounds.min_pt = center;
    ne_bounds.max_pt = node->bounds.max_pt;

    BoundingBox2d sw_bounds;
    sw_bounds.min_pt = node->bounds.min_pt;
    sw_bounds.max_pt = center;

    BoundingBox2d se_bounds;
    se_bounds.min_pt = Vec2d(center.x(), node->bounds.min_pt.y());
    se_bounds.max_pt = Vec2d(node->bounds.max_pt.x(), center.y());

    node->children[0] = make_unique_cpp11<Node>(nw_bounds); // NW
    node->children[1] = make_unique_cpp11<Node>(ne_bounds); // NE
    node->children[2] = make_unique_cpp11<Node>(sw_bounds); // SW
    node->children[3] = make_unique_cpp11<Node>(se_bounds); // SE

    node->divided = true;
}

bool QuadTree::intersects(const BoundingBox2d& box, const BezierCurve& curve) const {
    // Quick check using the curve's bounding box
    BoundingBox2d curve_box = curve.bbox();
    return box.intersects(curve_box);
}

BoundingBox2d QuadTree::getCurveBounds(const BezierCurve& curve) const {
    return curve.bbox();
}

void QuadTree::insert(const BezierCurve& curve, const std::string& id) {
    if (!root_->bounds.contains(getCurveBounds(curve).min_pt) ||
        !root_->bounds.contains(getCurveBounds(curve).max_pt)) {
        // Curve is outside the quadtree bounds
        return;
    }
    insertRecursive(root_.get(), curve, id, 0);
}

void QuadTree::insertRecursive(Node* node, const BezierCurve& curve, const std::string& id, size_t current_depth) {
    BoundingBox2d curve_bounds = getCurveBounds(curve);

    // If the node hasn't been subdivided yet and has room, or we've reached the depth limit
    if (!node->divided && node->curves.size() < capacity_ && current_depth < depth_limit_) {
        node->curves.emplace_back(curve, id);
        return;
    }

    // Subdivide if needed
    if (!node->divided && current_depth < depth_limit_) {
        subdivide(node);

        // Redistribute existing curves to children
        auto existing_curves = std::move(node->curves);
        for (const auto& curve_pair : existing_curves) {
            bool inserted = false;
            for (int i = 0; i < 4; i++) {
                if (node->children[i]->bounds.contains(getCurveBounds(curve_pair.first).min_pt) &&
                    node->children[i]->bounds.contains(getCurveBounds(curve_pair.first).max_pt)) {
                    insertRecursive(node->children[i].get(), curve_pair.first, curve_pair.second, current_depth + 1);
                    inserted = true;
                    break;
                }
            }
            if (!inserted) {
                // If it spans multiple quadrants, keep in current node
                node->curves.push_back(curve_pair);
            }
        }
    }

    // Insert the new curve
    bool inserted = false;
    if (current_depth < depth_limit_) {
        for (int i = 0; i < 4; i++) {
            if (node->children[i]->bounds.contains(curve_bounds.min_pt) &&
                node->children[i]->bounds.contains(curve_bounds.max_pt)) {
                insertRecursive(node->children[i].get(), curve, id, current_depth + 1);
                inserted = true;
                break;
            }
        }
    }

    if (!inserted) {
        // If it spans multiple quadrants or we've reached depth limit, keep in current node
        node->curves.emplace_back(curve, id);
    }
}

void QuadTree::queryRecursive(const Node* node, const BoundingBox2d& range,
                              std::vector<std::pair<BezierCurve, std::string>>& results) const {
    if (!node->bounds.intersects(range)) {
        return;
    }

    // Check curves in current node
    for (const auto& curve_pair : node->curves) {
        if (intersects(range, curve_pair.first)) {
            results.push_back(curve_pair);
        }
    }

    // If node is divided, check children
    if (node->divided) {
        for (int i = 0; i < 4; i++) {
            queryRecursive(node->children[i].get(), range, results);
        }
    }
}

std::vector<std::pair<BezierCurve, std::string>> QuadTree::queryRange(const BoundingBox2d& range) const {
    std::vector<std::pair<BezierCurve, std::string>> results;
    queryRecursive(root_.get(), range, results);
    return results;
}

void QuadTree::clear() {
    root_ = make_unique_cpp11<Node>(root_->bounds);
}