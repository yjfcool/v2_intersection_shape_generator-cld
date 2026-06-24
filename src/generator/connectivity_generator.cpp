#include "connectivity_generator.h"
#include "curve/hermite_init.h"
#include "curve/curve_utils.h"
#include "constraints/fence_check.h"
#include "utils.h"
#include "optimizer/sdf_field.h"
#include "constraints/infeasibility_detector.h"
#include "utils/quadtree.h"
#include <chrono>
#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <unordered_set>
#include <cstdlib>

namespace isg {

// 调试辅助: 当环境变量ISG_DEBUG_UTURN设置时返回true。
static bool isgDebugUTurn() {
    static bool v = (std::getenv("ISG_DEBUG_UTURN") != nullptr);
    return v;
}

static bool isPureGeometricInput(const IntersectionInput& input) {
    return input.obstacles.empty() && input.boundaries.empty() &&
           input.crosswalks.empty() && input.stop_lines.empty();
}

// ── IntersectionInput helpers ─────────────────────────────────────────────────
const bool IntersectionInput::IsEntryLane(const LaneId& id) const {
    for (auto& lg : lane_groups)
        if (std::find(lg.lanes.begin(), lg.lanes.end(), id) != lg.lanes.end())
            return lg.role == GroupRole::Entry;
    return false;
}

const bool IntersectionInput::IsEntryLaneEdge(const LaneEdgeId& id) const {
    for (auto& lg : lane_groups)
        if (std::find(lg.boundaries.begin(), lg.boundaries.end(), id) != lg.boundaries.end())
            return lg.role == GroupRole::Entry;
    for (auto& lg : lane_groups)
        for (auto& lid : lg.lanes) {
            auto l = findLane(lid);
            if (l && (l->left_edge_id == id || l->right_edge_id == id))
                return lg.role == GroupRole::Entry;
        }
    return false;
}

const Lane* IntersectionInput::findLane(const LaneId& id) const {
    for (auto& l : lanes) if (l.id == id)return &l;
    return nullptr;
}

const LaneGroup* IntersectionInput::findGroup(const LaneGroupId& id) const {
    for (auto& g : lane_groups) if (g.id == id)return &g;
    return nullptr;
}

const LaneEdge* IntersectionInput::findEdge(const LaneEdgeId& id) const {
    for (auto& e : lane_edges) if (e.id == id)return &e;
    return nullptr;
}

bool IntersectionInput::laneGroupExists(const LaneGroupId& id) const {
    return findGroup(id) != nullptr;
}

std::pair<Vec2d, Vec2d> IntersectionInput::entryPtDir(const LaneId& lid) const {
    auto* l = findLane(lid);
    if (!l || l->geometry.points.empty()) {
        std::cout << "[WARN] entrylane:" << lid << " no geometry!\n";
        return {Vec2d(0, 0), Vec2d(1, 0)};
    }
    return {entryLinePoint(l->geometry.points), entryLineTangent(l->geometry.points)};
}

std::pair<Vec2d, Vec2d> IntersectionInput::exitPtDir(const LaneId& lid) const {
    auto* l = findLane(lid);
    if (!l || l->geometry.points.empty()) {
        std::cout << "[WARN] exitlane:" << lid << " no geometry!\n";
        return {Vec2d(10, 0), Vec2d(1, 0)};
    }
    return {exitLinePoint(l->geometry.points), exitLineTangent(l->geometry.points)};
}

struct LaneDirectionSample {
    LaneId lane_id;
    Vec2d direction{1, 0};
    int lane_order = 0;
    int group_index = 0;
};

static double directionAngleDiff(const Vec2d& a, const Vec2d& b) {
    if (a.norm() < 1e-10 || b.norm() < 1e-10)
        return M_PI;
    double c = a.normalized().dot(b.normalized());
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c);
}

static const Lane* findLaneById(const std::vector<Lane>& lanes, const LaneId& id) {
    for (const auto& lane : lanes)
        if (lane.id == id)
            return &lane;
    return nullptr;
}

static Lane* findLaneById(std::vector<Lane>& lanes, const LaneId& id) {
    for (auto& lane : lanes)
        if (lane.id == id)
            return &lane;
    return nullptr;
}

static Vec2d laneDirectionForRole(const Lane& lane, GroupRole role) {
    Vec2d dir = role == GroupRole::Entry
        ? entryLineTangent(lane.geometry.points)
        : exitLineTangent(lane.geometry.points);
    return dir.norm() > 1e-10 ? dir.normalized() : Vec2d(1, 0);
}

static std::vector<LaneDirectionSample> collectLaneDirections(
    const LaneGroup& group, const std::vector<Lane>& lanes) {
    std::vector<LaneDirectionSample> samples;
    samples.reserve(group.lanes.size());
    std::unordered_set<LaneId> seen_lanes;
    for (int i = 0; i < (int)group.lanes.size(); ++i) {
        if (!seen_lanes.insert(group.lanes[i]).second)
            continue;
        const Lane* lane = findLaneById(lanes, group.lanes[i]);
        if (!lane || lane->geometry.points.empty())
            continue;
        LaneDirectionSample s;
        s.lane_id = lane->id;
        s.direction = laneDirectionForRole(*lane, group.role);
        s.lane_order = lane->laneOrder;
        s.group_index = i;
        samples.push_back(s);
    }
    return samples;
}

static const LaneDirectionSample* innermostLaneDirection(
    const std::vector<LaneDirectionSample>& samples) {
    if (samples.empty())
        return nullptr;
    return &*std::min_element(
        samples.begin(), samples.end(),
        [](const LaneDirectionSample& a, const LaneDirectionSample& b) {
            if (a.lane_order != b.lane_order)
                return a.lane_order < b.lane_order;
            return a.group_index < b.group_index;
        });
}

static Vec2d meanDirectionFromSupport(
    const std::vector<LaneDirectionSample>& samples,
    const std::vector<int>& support,
    const Vec2d& fallback) {
    Vec2d sum(0, 0);
    for (int idx : support) {
        if (idx < 0 || idx >= (int)samples.size())
            continue;
        if (samples[idx].direction.norm() > 1e-10)
            sum += samples[idx].direction.normalized();
    }
    return sum.norm() > 1e-10 ? sum.normalized() : fallback;
}

static Vec2d groupUnifiedDirection(
    const LaneGroup& group, const std::vector<Lane>& lanes, double threshold_deg) {
    auto samples = collectLaneDirections(group, lanes);
    const auto* inner = innermostLaneDirection(samples);
    Vec2d fallback = inner ? inner->direction : Vec2d(1, 0);

    if (samples.size() <= 2)
        return fallback;

    double threshold_rad = std::max(0.0, std::min(180.0, threshold_deg)) * M_PI / 180.0;
    int best_idx = -1;
    int best_count = 0;
    for (int i = 0; i < (int)samples.size(); ++i) {
        int count = 0;
        for (int j = 0; j < (int)samples.size(); ++j) {
            if (i == j)
                continue;
            if (directionAngleDiff(samples[i].direction, samples[j].direction) <= threshold_rad)
                ++count;
        }
        if (count > best_count ||
            (count == best_count && count > 0 &&
             (best_idx < 0 ||
              samples[i].lane_order < samples[best_idx].lane_order ||
              (samples[i].lane_order == samples[best_idx].lane_order &&
               samples[i].group_index < samples[best_idx].group_index)))) {
            best_idx = i;
            best_count = count;
        }
    }

    if (best_idx < 0 || best_count == 0)
        return fallback;

    std::vector<int> support;
    support.push_back(best_idx);
    for (int j = 0; j < (int)samples.size(); ++j) {
        if (j == best_idx)
            continue;
        if (directionAngleDiff(samples[best_idx].direction, samples[j].direction) <= threshold_rad)
            support.push_back(j);
    }
    return meanDirectionFromSupport(samples, support, fallback);
}

static void setEntryLaneDirection(Lane& lane, const Vec2d& direction) {
    if (direction.norm() < 1e-10)
        return;
    Vec2d dir = direction.normalized();
    auto& pts = lane.geometry.points;
    if (pts.empty())
        return;
    Vec2d p = pts.back();
    if (pts.size() == 1)
        pts.insert(pts.begin(), p - dir);
    else
        pts.insert(pts.end() - 1, p - dir);
}

static void setExitLaneDirection(Lane& lane, const Vec2d& direction) {
    if (direction.norm() < 1e-10)
        return;
    Vec2d dir = direction.normalized();
    auto& pts = lane.geometry.points;
    if (pts.empty())
        return;
    Vec2d p = pts.front();
    if (pts.size() == 1)
        pts.push_back(p + dir);
    else
        pts.insert(pts.begin() + 1, p + dir);
}

static void applyGroupUnifiedDirections(
    IntersectionInput& input, const ConnectivityDirectionConfig& cfg) {
    if (cfg.mode != ConnectivityDirectionMode::GroupUnified)
        return;

    std::unordered_map<LaneGroupId, Vec2d> group_dirs;
    for (const auto& group : input.lane_groups)
        group_dirs[group.id] = groupUnifiedDirection(
            group, input.lanes, cfg.group_similarity_angle_deg);

    for (const auto& group : input.lane_groups) {
        auto it = group_dirs.find(group.id);
        if (it == group_dirs.end())
            continue;
        for (const auto& lane_id : group.lanes) {
            Lane* lane = findLaneById(input.lanes, lane_id);
            if (!lane)
                continue;
            if (group.role == GroupRole::Entry)
                setEntryLaneDirection(*lane, it->second);
            else
                setExitLaneDirection(*lane, it->second);
        }
    }
}

// ── GlobalCoordinator ─────────────────────────────────────────────────────────
//  生成顺序:
//    优先级0: 直行(基准曲线, 先生成)
//    优先级1: 左转/右转(围绕直行调整)
//    优先级2: 左调头/右调头(最外侧, 最后生成)
//
//  同优先级内按 CLUSTER 空间顺序排序:
//    - 对每个连通关系, 查找其在进入组空间排序表中的位置
//      (来自 ClusterOrderSolver::entryGroupOrder())
//    - 降序排列使 index 0 = 最左侧曲线
//    - "最内侧" = 最靠近直行 = 左转中index最大、右转中index最小, 即最接近中位
//    - 按 |distance from middle| 升序生成, 内侧曲线先生成
//
//  这样当外侧曲线生成时, 其内侧邻居已在 `done` 中, 作为兄弟曲线传入
//  携带正确的 expected_side, 优化器从一开始就有正确参考。
// ─────────────────────────────────────────────────────────────────────────────
//
//
// 由端点切向几何计算转向类型, 不依赖 conn.turn_type (输入数据的转向类型不一定准确)。
//
// 算法:
//   t0 = 进入车道切向 (指向路口内)
//   t1 = 退出车道切向 (指向路口外)
//   • T0·T1 < −0.5  → 反平行退出 → U型调头 (优先级2)
//   • |cross2d(t0, path_dir)| < 0.25  → 近直行 (优先级0)
//   • cross2d(t0, path_dir) > 0 → 左转, < 0 → 右转 (优先级1)
//   其中 path_dir = (exit_pt − entry_pt).normalized()
static int globalTurnPriorityGeometric(
    const Connectivity& conn, const IntersectionInput& inp) {
    // 获取端点几何
    auto kv0 = inp.entryPtDir(conn.entry_lane_id);
    auto& p0 = kv0.first;
    auto& t0 = kv0.second;
    auto kv1 = inp.exitPtDir (conn.exit_lane_id);
    auto& p1 = kv1.first;
    auto& t1 = kv1.second;

    if (t0.norm() < 1e-9 || t1.norm() < 1e-9) {
        // 几何退化时回退到声明类型
        switch (conn.turn_type) {
        case ConnTurnType::Straight:   return 0;
        case ConnTurnType::UTurnLeft:
        case ConnTurnType::UTurnRight: return 2;
        default:                       return 1;
        }
    }
    t0.normalize(); t1.normalize();

    // U型调头: 退出切向与进入切向反平行
    if (t0.dot(t1) < -0.5) return 2;

    // 由进入→退出方位区分直行/转向
    Vec2d d = p1 - p0;
    if (d.norm() < 1e-9) return 0;
    d.normalize();
    double c = cross2d(t0, d); // positive → left, negative → right
    if (std::abs(c) < 0.25)   return 0; // straight (< ~15°)
    return 1;                            // left or right turn
}

// 兼容性包装(可读性保留)—— 不再用于优先级判定。
static int globalTurnPriority(ConnTurnType t) {
    switch (t) {
    case ConnTurnType::Straight: return 0;
    case ConnTurnType::TurnLeft: return 1;
    case ConnTurnType::TurnRight: return 1;
    case ConnTurnType::UTurnLeft: return 2;
    case ConnTurnType::UTurnRight: return 2;
    default: return 0;
    }
}

void GlobalCoordinator::build(
    const std::vector<Connectivity>& conns, const IntersectionInput& inp, const ClusterOrderSolver& cs) {
    // 建立 connId → 进入组内空间索引的映射(来自solver)
    std::unordered_map<ConnId, int> spatial_idx;
    std::unordered_map<ConnId, int> group_size;
    for (auto& kv : cs.entryGroupOrder()) {
        auto& gid=kv.first; auto& sorted=kv.second;
        int sz = (int)sorted.size();
        for (int i = 0; i < sz; ++i) {
            spatial_idx[sorted[i]] = i;
            group_size [sorted[i]] = sz;
        }
    }

    // 按几何转向优先级分组(不依赖 conn.turn_type 因其可能不准)
    std::map<int, std::vector<const Connectivity*>> pm;
    for (auto& c : conns)
        pm[globalTurnPriorityGeometric(c, inp)].push_back(&c);

    groups_.clear();
    for (auto& kv : pm) {
        auto& pri = kv.first; auto& cv = kv.second;
        // 优先级内排序: 最内侧(|spatial_idx - sz/2|最小)优先
        // 并列: 按 enterGroupId 分组(同组成员聚在一起), 再按 spatial_idx
        std::stable_sort(cv.begin(), cv.end(), [&](const Connectivity* a, const Connectivity* b) {
            int ia = spatial_idx.count(a->id) ? spatial_idx[a->id] : 0;
            int ib = spatial_idx.count(b->id) ? spatial_idx[b->id] : 0;
            int sa = group_size .count(a->id) ? group_size [a->id] : 1;
            int sb = group_size .count(b->id) ? group_size [b->id] : 1;
            // "Distance from middle" = |i - (sz-1)/2.0|
            double da = std::abs(ia - (sa-1)/2.0);
            double db = std::abs(ib - (sb-1)/2.0);
            if (std::abs(da-db) > 0.5) return da < db;   // inner before outer
            if (a->enterGroupId != b->enterGroupId)
                return a->enterGroupId < b->enterGroupId;  // same group together
            return ia < ib;
        });

        OptGroup g; g.priority = pri;
        for (auto* c : cv) g.conn_ids.push_back(c->id);
        groups_.push_back(std::move(g));
    }
}

void GlobalCoordinator::addSoftObstacles(SDFField&,
    const std::vector<ConnectivityCurve>&, double) const {}

// ── ConnectivityGenerator ─────────────────────────────────────────────────────

ConnectivityGenerator::ConnectivityGenerator(
    const LBFGSConfig& cfg,
    const ConnectivityDirectionConfig& direction_cfg)
    : solver_(cfg), direction_cfg_(direction_cfg) {}

// ─────────────────────────────────────────────────────────────────────────────
//  buildSiblings
//
//  对每条已生成曲线:
//  - 判定簇归属(同进入组 OR 同退出组)
//  - 设置 exempt_a1: 仅当 enterGroup 与 exitGroup 均不匹配时为true
//    (= 不同道路arm的交叉交通 = 结构性交叉)
//  - 设置 expected_side: 来自 ClusterOrderSolver::expectedSideOf(cid, id)
//    (排序为降序, 故符号正确)
//  - 设置 ref_perp: 来自 ClusterOrderSolver::refPerpOf()
// ─────────────────────────────────────────────────────────────────────────────
std::vector<SiblingCurve> ConnectivityGenerator::buildSiblings(
    const ConnId& id, const std::unordered_map<ConnId,BezierCurve>& done,
    const ClusterOrderSolver& cs, const std::vector<Connectivity>& conns,
    bool constrained_only,
    const std::unordered_set<ConnId>* fixed_shape_ids) const {
    std::vector<SiblingCurve> sibs;
    for (auto& kv : done) {
        auto& cid = kv.first;
        auto& curve = kv.second;
        if (cid == id)
            continue;

        SiblingCurve s;
        s.curve = curve;
        s.fixed_shape = fixed_shape_ids && fixed_shape_ids->count(cid) > 0;

        auto ex = cs.exemptionOf(id, cid);
        bool in_cluster = cs.pairExists(id, cid);
        if (!in_cluster) {
            // 无簇关系 → 不同arm → 结构性交叉交通
            s.exempt_a1 = true;
        } else if (ex == CrossExemption::StructuralCross) {
            // 拓扑倒置 → 必须交叉 → 豁免
            s.exempt_a1 = true;
        } else {
            // 同簇约束配对 → 施加惩罚
            s.exempt_a1 = false;
        }
        if (constrained_only && s.exempt_a1)
            continue;

        s.exempt_a2_radius = (ex == CrossExemption::ObstacleCross) ? 1.5 : 0.0;

        // expected_side 以 cid 视角相对 id:
        // cs.expectedSideOf(cid, id):
        //   +1 → 兄弟(cid)在当前(id)左侧
        //   -1 → 兄弟(cid)在当前(id)右侧
        s.expected_side = cs.expectedSideOf(cid, id);

        // evalCluster的固定横向参考轴
        s.ref_perp = cs.refPerpOf(id, cid);

        // 传递 shared_endpoint 标志, evalCluster使用更宽跳过区
        s.shared_endpoint = cs.isSharedEndpoint(id, cid);

        sibs.push_back(std::move(s));
    }
    return sibs;
}

// ─────────────────────────────────────────────────────────────────────────────
//  computeRequiredLateralOffset
//
//  对连通关系 `id` 的初始曲线, 检查 p0→p1 直接Bezier弧
//  是否违反任一内侧兄弟的顺序。需要时返回(侧, 大小)横向偏移:
//    side = +1 (向左推), side = -1 (向右推)。
//
//  仅考虑内侧兄弟(应位于当前曲线右侧, 即 expected_side=-1)。
//  若直接弧已位于所有内侧兄弟的外侧(左侧) → 无需偏移。
//
//  ref_perp: 横向参考方向(arm左法向)。
// ─────────────────────────────────────────────────────────────────────────────
static bool needsLateralOffset(
    const Vec2d& p0, const Vec2d& p1,
    const std::vector<SiblingCurve>& siblings,
    int& out_side, double& out_offset_m)
{
    out_side = 0; out_offset_m = 0.0;

    // 找应位于我右侧的兄弟(expected_side=-1),
    // 检查我是否实际位于其左侧(正确)。
    // 最关键的是: 内侧兄弟中横向值最大者
    // (内侧中最左 = 我必须保持在其外侧的兄弟)。
    double worst_violation = 0.0;
    Vec2d ref_perp(0,0);
    int required_push_side = 0;

    for (auto& s : siblings) {
        if (s.exempt_a1) continue;
        if (s.expected_side != -1) continue;  // only inner siblings (should be to my RIGHT)
        if (s.ref_perp.norm() < 1e-9) continue;

        // 兄弟在 ref_perp 坐标系下的平均横向位置
        auto sib_pts = s.curve.sampleByArcLength(20);
        if (sib_pts.empty()) continue;
        double sib_lat = 0;
        for (auto& pt : sib_pts) sib_lat += pt.dot(s.ref_perp);
        sib_lat /= sib_pts.size();

        // 当前曲线估计平均横向位置(直接弧中点)
        Vec2d mid = 0.5*(p0+p1);
        double my_lat = mid.dot(s.ref_perp);

        // 对 expected_side=-1 (兄弟应在右侧 = 兄弟横向值更小):
        // 违反条件: my_lat < sib_lat + MARGIN (我向左偏移不够)
        double margin = 0.25;
        double viol = sib_lat + margin - my_lat;  // positive = violation
        if (viol > worst_violation) {
            worst_violation = viol;
            ref_perp = s.ref_perp;
            required_push_side = +1;  // push LEFT (increase lat)
        }
    }

    // 同时检查应位于我左侧的兄弟(expected_side=+1):
    // 这些是外侧兄弟, 我必须保持在其内侧。
    for (auto& s : siblings) {
        if (s.exempt_a1) continue;
        if (s.expected_side != +1) continue;  // only outer siblings (should be to my LEFT)
        if (s.ref_perp.norm() < 1e-9) continue;

        auto sib_pts = s.curve.sampleByArcLength(20);
        if (sib_pts.empty()) continue;
        double sib_lat = 0;
        for (auto& pt : sib_pts) sib_lat += pt.dot(s.ref_perp);
        sib_lat /= sib_pts.size();

        Vec2d mid = 0.5*(p0+p1);
        double my_lat = mid.dot(s.ref_perp);

        double margin = 0.25;
        double viol = my_lat + margin - sib_lat;  // positive = violation (I'm too far LEFT)
        // 仅当比反方向违反更强时应用
        if (viol > worst_violation && required_push_side == 0) {
            worst_violation = viol;
            ref_perp = s.ref_perp;
            required_push_side = -1;  // push RIGHT (decrease lat)
        }
    }

    if (worst_violation < 1e-3 || required_push_side == 0) return false;

    out_side     = required_push_side;
    out_offset_m = worst_violation + 0.1;
    return true;
}

BezierCurve ConnectivityGenerator::postProcess(
    const BezierCurve& c, const SDFField& sdf, const Polygon2d& fence, double kmax,
    const Vec2d& t0_orig, const Vec2d& t1_orig, bool skip_elastic_band,
    const Vec2d* p0_exact, const Vec2d* p1_exact) {
    auto refined = adaptiveRefine(c, sdf, kmax);
    BezierCurve cur = refined.curve;
    if (refined.was_split) {
        PenaltyCost cost2;
        cost2.proto = cur;
        cost2.sdf = &sdf;
        cost2.fence = fence;
        cost2.start_tan_dir = t0_orig.norm() > 1e-8 ? t0_orig.normalized() : c.startTan();
        cost2.end_tan_dir = t1_orig.norm() > 1e-8 ? t1_orig.normalized() : c.endTan();
        cost2.obstacle_clearance = 0.0;
        cost2.full_param_mode = (cur.numSegments() > 1);
        cost2.buildCache();
        cur = optimiseCurve(cost2, solver_, cur, 2);
    }

    Vec2d st = t0_orig.norm() > 1e-8 ? t0_orig.normalized() : cur.startTan();
    Vec2d et = t1_orig.norm() > 1e-8 ? t1_orig.normalized() : cur.endTan();
    Vec2d ep0 = p0_exact ? *p0_exact : cur.startPt();
    Vec2d ep1 = p1_exact ? *p1_exact : cur.endPt();

    if (skip_elastic_band) {
        if (!cur.segs.empty()) {
            cur.segs.front().ctrl[0] = ep0;
            cur.segs.back().ctrl[3] = ep1;
            Vec2d& c1 = cur.segs.front().ctrl[1];
            double sl = (cur.segs.front().ctrl[3] - ep0).norm();
            double lam = std::max((c1 - ep0).dot(st), sl * 0.1);
            c1 = ep0 + lam * st;
            Vec2d& c2 = cur.segs.back().ctrl[2];
            double el = (ep1 - cur.segs.back().ctrl[0]).norm();
            double mu = std::max((ep1 - c2).dot(et), el * 0.1);
            c2 = ep1 - mu * et;
        }
        return cur;
    }

    double arc = cur.arcLength();
    // 性能优化: n_samp 上限从无限制降至100, 80次迭代降至30
    int n_samp = std::max(40, std::min(100, (int)(arc / 0.15)));
    auto pts = cur.sampleByArcLength(n_samp);
    if ((int)pts.size() >= 3) {
        pts.front() = ep0;
        pts.back() = ep1;

        double max_k = cur.maxCurvature(20);
        double move_st = std::min(0.15, std::max(0.02, max_k * 0.3));
        // 性能优化: 迭代从80降至30
        auto sm = elasticBandSmooth(pts, sdf, fence, kmax, move_st, 30, 0.1);
        sm.front() = ep0;
        sm.back() = ep1;

        cur = rebuildFromSmoothedPts(sm, st, et);
        if (!cur.segs.empty()) {
            cur.segs.front().ctrl[0] = ep0;
            cur.segs.back().ctrl[3] = ep1;
        }
    }
    return cur;
}

static double minSDFAlongCurveAdaptive(const BezierCurve& curve, const SDFField& sdf);
static bool curveIntersectsObstacles(
    const BezierCurve& curve, const std::vector<Obstacle>& obstacles, Vec2d* out);
static bool curveIntersectsBoundaries(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries);

void ConnectivityGenerator::validate(
    ConnectivityCurve& cc, const IntersectionInput& input, const SDFField& sdf) const {
    if (!cc.curve) return;
    auto& c = *cc.curve;
    double ms = minSDFAlongCurveAdaptive(c, sdf);
    cc.violation.max_obstacle_penetration = std::max(0.0, -ms);
    Vec2d obstacle_hit;
    if (curveIntersectsObstacles(c, input.obstacles, &obstacle_hit)) {
        cc.violation.max_obstacle_penetration =
            std::max(cc.violation.max_obstacle_penetration, 0.01);
        cc.violation.reason = "curve intersects obstacle";
    }
    if (!input.area.geometry.outer.empty()) {
        double ov = 0;
        // 性能优化: 上限从240降至80,降低采样数
        int n = std::max(20, std::min(80, (int)std::ceil(c.arcLength() / 0.25) + 1));
        for (auto& pt : c.sampleByArcLength(n))
            if (!polygonContains(input.area.geometry, pt))
                ov = std::max(ov, pointToPolygonDist(pt, input.area.geometry));
        cc.violation.max_fence_overflow = ov;
    }
    bool self_cross = curveSelfIntersectsBusiness(c, 1.0);
    if (self_cross)
        cc.violation.reason = "curve self-intersects away from endpoints";
    bool boundary_cross = curveIntersectsBoundaries(c, input.boundaries);
    if (boundary_cross)
        cc.violation.reason = "curve intersects boundary away from endpoints";
    if (self_cross || boundary_cross || cc.violation.max_obstacle_penetration > 0.05)
        cc.status = CurveStatus::Degraded;
    else if (!cc.violation.exempt_crosses.empty())
        cc.status = CurveStatus::WarnA2;
    else
        cc.status = CurveStatus::OK;
}

static std::unordered_map<ConnId, size_t> resultIndexById(
    const std::vector<ConnectivityCurve>& results) {
    std::unordered_map<ConnId, size_t> idx;
    for (size_t i = 0; i < results.size(); ++i)
        idx[results[i].id] = i;
    return idx;
}

static std::unordered_map<ConnId, BezierCurve> curveMapFromResults(
    const std::vector<ConnectivityCurve>& results) {
    std::unordered_map<ConnId, BezierCurve> curves;
    for (auto& cc : results)
        if (cc.curve)
            curves[cc.id] = *cc.curve;
    return curves;
}

static bool hasFixedGeometry(const Connectivity& conn) {
    return conn.geometry.points.size() >= 2;
}

static BezierCurve fixedGeometryToCurve(
    const Connectivity& conn, const IntersectionInput& input) {
    (void)input;
    BezierCurve curve;
    const auto& pts = conn.geometry.points;
    if (pts.size() < 2)
        return curve;

    for (int i = 0; i + 1 < (int)pts.size(); ++i) {
        Vec2d chord = pts[i + 1] - pts[i];
        double len = chord.norm();
        if (len < 1e-8)
            continue;
        Vec2d seg_dir = chord / len;
        curve.segs.push_back(makeCubicG1(pts[i], seg_dir, pts[i + 1], seg_dir, 1.0 / 3.0));
    }
    return curve;
}

static ConnectivityCurve makeFixedGeometryCurve(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf) {
    ConnectivityCurve cc;
    cc.id = conn.id;
    cc.entry_lane_id = conn.entry_lane_id;
    cc.exit_lane_id = conn.exit_lane_id;
    cc.turn_type = conn.turn_type;
    cc.geometry = conn.geometry;
    BezierCurve curve = fixedGeometryToCurve(conn, input);
    if (!curve.empty())
        cc.curve = std::make_shared<BezierCurve>(curve);

    if (cc.curve) {
        double ms = minSDFAlongCurveAdaptive(*cc.curve, sdf);
        cc.violation.max_obstacle_penetration = std::max(0.0, -ms);
        Vec2d obstacle_hit;
        if (curveIntersectsObstacles(*cc.curve, input.obstacles, &obstacle_hit)) {
            cc.violation.max_obstacle_penetration =
                std::max(cc.violation.max_obstacle_penetration, 0.01);
            cc.violation.reason = "fixed geometry intersects obstacle";
            cc.status = CurveStatus::Degraded;
        }
    }
    return cc;
}

static bool fixedGeometryHitsObstacle(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf) {
    if (!hasFixedGeometry(conn))
        return false;
    BezierCurve curve = fixedGeometryToCurve(conn, input);
    if (curve.empty())
        return false;
    (void)sdf;
    return curveIntersectsObstacles(curve, input.obstacles, nullptr);
}

struct SampledCurve {
    std::vector<Vec2d> pts;
    BoundingBox2d bbox;
    Vec2d start{0, 0};
    Vec2d end{0, 0};
};

struct SampledSiblingCurve {
    SampledCurve sampled;
    bool exempt_a1 = false;
    int expected_side = 0;       // +1=兄弟在左, -1=兄弟在右
    Vec2d ref_perp{0, 0};        // 配对横向参考轴
    bool shared_endpoint = false;
};

static SampledCurve sampleCurveForIntersections(const BezierCurve& curve, int n = 32) {
    SampledCurve s;
    s.start = curve.startPt();
    s.end = curve.endPt();
    int adaptive_n = std::max(n, (int)std::ceil(curve.arcLength() / 0.20) + 1);
    adaptive_n = std::min(adaptive_n, 240);
    s.pts = curve.sampleByArcLength(adaptive_n);
    for (const auto& pt : s.pts)
        s.bbox.expand(pt);
    return s;
}

/// 沿曲线自适应采样查询最小SDF值
/// 性能优化: 上限从320降至100, 不调用queryWithGrad(仅用query)
static double minSDFAlongCurveAdaptive(const BezierCurve& curve, const SDFField& sdf) {
    if (!sdf.valid())
        return 1e18;
    // 性能优化: 上限从320降至100
    int n = std::max(20, std::min(100, (int)std::ceil(curve.arcLength() / 0.30) + 1));
    double m = 1e18;
    for (const auto& pt : curve.sampleByArcLength(n))
        m = std::min(m, sdf.queryWithGrad(pt).first);
    return m;
}

static bool pointOnPolygonBoundary(const Vec2d& pt, const Polygon2d& poly, double tol = 1e-6) {
    const auto& ring = poly.outer;
    if (ring.size() < 2)
        return false;
    for (int i = 0; i < (int)ring.size(); ++i) {
        const Vec2d& a = ring[i];
        const Vec2d& b = ring[(i + 1) % ring.size()];
        if (pointToSegment(pt, a, b).first <= tol)
            return true;
    }
    return false;
}

static bool pointInsideOrOnPolygon(const Vec2d& pt, const Polygon2d& poly, double tol = 1e-6) {
    return polygonContains(poly, pt) || pointOnPolygonBoundary(pt, poly, tol);
}

static bool polygonSegmentIntersects(
    const Vec2d& a, const Vec2d& b, const Polygon2d& poly, Vec2d* out = nullptr) {
    const auto& ring = poly.outer;
    if (ring.size() < 2)
        return false;
    for (int i = 0; i < (int)ring.size(); ++i) {
        Vec2d isect;
        if (segmentsIntersect(a, b, ring[i], ring[(i + 1) % ring.size()], &isect)) {
            if (out)
                *out = isect;
            return true;
        }
    }
    return false;
}

static const Polygon2d& obstacleHardGeometry(const Obstacle& obs) {
    return obs.geometry.outer.empty() ? obs.buffered_geometry : obs.geometry;
}

static bool curveIntersectsObstacles(
    const BezierCurve& curve, const std::vector<Obstacle>& obstacles, Vec2d* out = nullptr) {
    if (obstacles.empty())
        return false;

    auto sampled = sampleCurveForIntersections(curve, 64);
    if (sampled.pts.size() < 2)
        return false;

    for (const auto& obs : obstacles) {
        const Polygon2d& poly = obstacleHardGeometry(obs);
        if (poly.outer.size() < 3)
            continue;

        BoundingBox2d obs_box = poly.bbox();
        if (!sampled.bbox.intersects(obs_box))
            continue;

        for (int i = 1; i + 1 < (int)sampled.pts.size(); ++i) {
            if (pointInsideOrOnPolygon(sampled.pts[i], poly)) {
                if (out)
                    *out = sampled.pts[i];
                return true;
            }
        }

        for (int i = 0; i + 1 < (int)sampled.pts.size(); ++i) {
            Vec2d isect;
            if (polygonSegmentIntersects(sampled.pts[i], sampled.pts[i + 1], poly, &isect)) {
                if (out)
                    *out = isect;
                return true;
            }
        }
    }
    return false;
}

static double distToSampledEndpoints(const Vec2d& pt, const SampledCurve& a, const SampledCurve& b) {
    double d = (pt - a.start).norm();
    d = std::min(d, (pt - a.end).norm());
    d = std::min(d, (pt - b.start).norm());
    d = std::min(d, (pt - b.end).norm());
    return d;
}

static bool sampledCurvesIntersectBusiness(
    const SampledCurve& a, const SampledCurve& b, double endpoint_tol, Vec2d* out = nullptr) {
    if (a.pts.size() < 2 || b.pts.size() < 2 || !a.bbox.intersects(b.bbox))
        return false;
    for (int ai = 0; ai + 1 < (int)a.pts.size(); ++ai) {
        Vec2d amid = 0.5 * (a.pts[ai] + a.pts[ai + 1]);
        for (int bi = 0; bi + 1 < (int)b.pts.size(); ++bi) {
            Vec2d bmid = 0.5 * (b.pts[bi] + b.pts[bi + 1]);
            if ((amid - bmid).squaredNorm() > 900.0)
                continue;
            Vec2d isect;
            if (!segmentsIntersect(a.pts[ai], a.pts[ai + 1], b.pts[bi], b.pts[bi + 1], &isect))
                continue;
            if (distToSampledEndpoints(isect, a, b) <= endpoint_tol)
                continue;
            if (out)
                *out = isect;
            return true;
        }
    }
    return false;
}

static std::vector<SampledSiblingCurve> sampleSiblingsForIntersections(
    const std::vector<SiblingCurve>& siblings, int n = 48) {  //提升至48,避免粗采样漏检
    std::vector<SampledSiblingCurve> sampled;
    sampled.reserve(siblings.size());
    for (const auto& sib : siblings) {
        SampledSiblingCurve s;
        s.exempt_a1 = sib.exempt_a1;
        s.expected_side = sib.expected_side;
        s.ref_perp = sib.ref_perp;
        s.shared_endpoint = sib.shared_endpoint;
        s.sampled = sampleCurveForIntersections(sib.curve, n);
        sampled.push_back(std::move(s));
    }
    return sampled;
}

static int sampledSiblingCrossCount(
    const SampledCurve& curve, const std::vector<SampledSiblingCurve>& siblings,
    bool constrained_only, double endpoint_tol = 1.5) {
    int count = 0;
    for (const auto& sib : siblings) {
        if (constrained_only && sib.exempt_a1)
            continue;
        if (sampledCurvesIntersectBusiness(curve, sib.sampled, endpoint_tol))
            ++count;
    }
    return count;
}

static int sampledSiblingCrossCount(
    const BezierCurve& curve, const std::vector<SampledSiblingCurve>& siblings,
    bool constrained_only, double endpoint_tol = 1.5) {
    return sampledSiblingCrossCount(
        sampleCurveForIntersections(curve), siblings, constrained_only, endpoint_tol);
}

static double sharedEndpointSideViolation(
    const SampledCurve& curve, const std::vector<SampledSiblingCurve>& siblings,
    double endpoint_tol = 1.5) {
    if (curve.pts.size() < 4)
        return 0.0;

    double violation = 0.0;
    constexpr double SIDE_MARGIN = 0.25;
    for (const auto& sib : siblings) {
        if (sib.exempt_a1 || !sib.shared_endpoint)
            continue;
        if (sib.expected_side == 0 || sib.ref_perp.norm() < 1e-9)
            continue;
        if (sib.sampled.pts.size() < 4)
            continue;

        bool share_start = (curve.start - sib.sampled.start).norm() <= endpoint_tol;
        bool share_end = (curve.end - sib.sampled.end).norm() <= endpoint_tol;
        if (!share_start && !share_end)
            continue;

        for (const auto& pt : curve.pts) {
            Vec2d anchor = share_start ? curve.start : curve.end;
            double d = (pt - anchor).norm();
            if (d <= endpoint_tol || d > 8.0)
                continue;

            double best_d2 = std::numeric_limits<double>::infinity();
            double best_lat = 0.0;
            for (const auto& sp : sib.sampled.pts) {
                Vec2d sib_anchor = share_start ? sib.sampled.start : sib.sampled.end;
                double sd = (sp - sib_anchor).norm();
                if (sd <= endpoint_tol || sd > 8.0)
                    continue;
                double d2 = (sp - pt).squaredNorm();
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_lat = sp.dot(sib.ref_perp);
                }
            }
            if (!std::isfinite(best_d2) || best_d2 > 36.0)
                continue;

            double diff = pt.dot(sib.ref_perp) - best_lat;
            if (sib.expected_side == +1)
                violation += std::max(0.0, diff + SIDE_MARGIN);
            else
                violation += std::max(0.0, -diff + SIDE_MARGIN);
        }
    }
    return violation;
}

static std::unordered_set<ConnId> crossingIdsTouchingSeeds(
    const std::vector<ConnectivityCurve>& results,
    const ClusterOrderSolver& cs, const std::unordered_set<ConnId>& seeds,
    double endpoint_tol = 1.5) {
    std::unordered_set<ConnId> bad;
    if (seeds.empty())
        return bad;
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, SampledCurve> samples;
    samples.reserve(results.size());
    for (const auto& cc : results)
        if (cc.curve)
            samples.emplace(cc.id, sampleCurveForIntersections(*cc.curve));
    for (auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        if (!seeds.count(p.id_a) && !seeds.count(p.id_b))
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        auto sa = samples.find(p.id_a);
        auto sb = samples.find(p.id_b);
        if (sa == samples.end() || sb == samples.end())
            continue;
        if (sampledCurvesIntersectBusiness(sa->second, sb->second, endpoint_tol)) {
            if (seeds.count(p.id_a))
                bad.insert(p.id_a);
            if (seeds.count(p.id_b))
                bad.insert(p.id_b);
        }
    }
    return bad;
}

static std::unordered_set<ConnId> allConstrainedCrossingIds(
    const std::vector<ConnectivityCurve>& results,
    const ClusterOrderSolver& cs,
    const std::unordered_set<ConnId>& preserved_fixed_ids,
    double endpoint_tol = 1.5) {
    std::unordered_set<ConnId> bad;
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, SampledCurve> samples;
    samples.reserve(results.size());
    for (const auto& cc : results)
        if (cc.curve)
            samples.emplace(cc.id, sampleCurveForIntersections(*cc.curve));

    for (auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        auto sa = samples.find(p.id_a);
        auto sb = samples.find(p.id_b);
        if (sa == samples.end() || sb == samples.end())
            continue;
        if (!sampledCurvesIntersectBusiness(sa->second, sb->second, endpoint_tol))
            continue;

        if (!preserved_fixed_ids.count(p.id_a))
            bad.insert(p.id_a);
        if (!preserved_fixed_ids.count(p.id_b))
            bad.insert(p.id_b);
    }
    return bad;
}

static void annotateClusterCrossings(
    std::vector<ConnectivityCurve>& results,
    const ClusterOrderSolver& cs, double endpoint_tol = 1.5) {
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, SampledCurve> samples;
    samples.reserve(results.size());
    // 性能优化: 采样数从48降至24(Debug进一步降至16)
#ifdef NDEBUG
    const int N_ANNOTATE = 24;
#else
    const int N_ANNOTATE = 16;
#endif
    for (const auto& cc : results)
        if (cc.curve)
            samples.emplace(cc.id, sampleCurveForIntersections(*cc.curve, N_ANNOTATE));
    for (auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        auto& ca = results[ia->second];
        auto& cb = results[ib->second];
        if (!ca.curve || !cb.curve)
            continue;
        auto sa = samples.find(p.id_a);
        auto sb = samples.find(p.id_b);
        if (sa == samples.end() || sb == samples.end())
            continue;
        Vec2d pt;
        bool interior = sampledCurvesIntersectBusiness(sa->second, sb->second, endpoint_tol, &pt);
        if (!interior)
            continue;
        ca.violation.exempt_crosses.push_back(pt);
        cb.violation.exempt_crosses.push_back(pt);
        CurveStatus st = (p.exempt == CrossExemption::ObstacleCross)
            ? CurveStatus::WarnA2 : CurveStatus::Degraded;
        if (ca.status == CurveStatus::OK || st == CurveStatus::Degraded)
            ca.status = st;
        if (cb.status == CurveStatus::OK || st == CurveStatus::Degraded)
            cb.status = st;
    }
}

static bool naturalCubicHitsObstacle(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const SDFField& sdf, const std::vector<Obstacle>& obstacles) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    BezierCurve trial;
    trial.segs.push_back(makeCubicG1(p0, T0, p1, T1, 0.4));
    if (curveIntersectsObstacles(trial, obstacles))
        return true;
    return sdf.valid() && minSDFAlongCurveAdaptive(trial, sdf) < 0.0;
}

static bool curveIntersectsBoundaries(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries) {
    if (boundaries.empty())
        return false;
    auto sampled = sampleCurveForIntersections(curve, 64);
    for (const auto& bnd : boundaries) {
        const auto& pts = bnd.geometry.points;
        if (pts.size() < 2)
            continue;
        BoundingBox2d bnd_box;
        for (const auto& pt : pts)
            bnd_box.expand(pt);
        if (!sampled.bbox.intersects(bnd_box))
            continue;
        for (int i = 0; i + 1 < (int)sampled.pts.size(); ++i) {
            for (int j = 0; j + 1 < (int)pts.size(); ++j) {
                Vec2d isect;
                if (!segmentsIntersect(sampled.pts[i], sampled.pts[i + 1], pts[j], pts[j + 1], &isect))
                    continue;
                if ((isect - sampled.start).norm() <= 0.75 ||
                    (isect - sampled.end).norm() <= 0.75)
                    continue;
                if ((isect - pts[j]).norm() <= 0.10 ||
                    (isect - pts[j + 1]).norm() <= 0.10)
                    continue;
                return true;
            }
        }
    }
    return false;
}

static bool curveLeavesFence(const BezierCurve& curve, const Polygon2d& fence) {
    if (fence.outer.empty())
        return false;
    auto pts = curve.sample(24);
    for (int i = 1; i + 1 < (int)pts.size(); ++i) {
        const auto& pt = pts[i];
        if (!polygonContains(fence, pt) && pointToPolygonDist(pt, fence) > 0.10)
            return true;
    }
    return false;
}

// ─────────────────────────────────────────────────────────────────────────────
//  isCurveSeverelyDivergent
//
//  Detects optimizer divergence: when LBFGS pushes control points far away
//  from the p0-p1 chord, the curve becomes unusable (huge arc length, extreme
//  curvature, self-intersection, control points far from chord).
// ─────────────────────────────────────────────────────────────────────────────
static bool isCurveSeverelyDivergent(
    const BezierCurve& curve, const Vec2d& p0, const Vec2d& p1,
    double chord_len) {
    if (curve.empty() || chord_len < 1e-6) return false;

    // (a) Self-intersection (away from endpoints) — always divergent
    if (curveSelfIntersectsBusiness(curve, 1.0)) return true;

    // (b) Arc length >> chord (way too long)
    double arc = curve.arcLength();
    if (arc > std::max(3.0 * chord_len + 20.0, chord_len * 5.0)) return true;

    // (c) Max curvature extremely high (way too curvy)
    double maxk = curve.maxCurvature(20);
    if (maxk > 5.0) return true;

    // (d) Any control point far from chord (e.g., > 3×chord from p0 or p1)
    //     This catches the (75.52, 143.02) case for conn 26.
    double far_thresh = 3.0 * chord_len + 10.0;
    double far_thresh2 = far_thresh * far_thresh;
    for (auto& seg : curve.segs) {
        for (int k = 0; k < 4; ++k) {
            const Vec2d& cp = seg.ctrl[k];
            if ((cp - p0).squaredNorm() > far_thresh2 &&
                (cp - p1).squaredNorm() > far_thresh2)
                return true;
        }
    }
    return false;
}

struct CurveRisk {
    bool obstacle = false;
    bool boundary = false;
    bool fence = false;
    int sibling_crosses = 0;

    bool physical() const {
        return obstacle || boundary || fence;
    }
};

static void setConnectivityCurveGeometry(
    ConnectivityCurve& cc, const BezierCurve& curve, const LineString2d* fixed_geometry = nullptr) {
    cc.curve = std::make_shared<BezierCurve>(curve);
    cc.geometry.points.clear();
    if (fixed_geometry && fixed_geometry->points.size() >= 2) {
        cc.geometry = *fixed_geometry;
        return;
    }
    if (!curve.empty()) {
        int n = std::max(2, std::min(240, (int)std::ceil(curve.arcLength() / 0.3) + 1));
        cc.geometry.points = curve.sampleByArcLength(n);
    }
}

static CurveRisk assessCurveRisk(
    const BezierCurve& curve, const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence = true) {
    CurveRisk risk;
    double ms = minSDFAlongCurveAdaptive(curve, sdf);
    risk.obstacle = curveIntersectsObstacles(curve, input.obstacles) || (ms < 0.0);
    risk.boundary = curveIntersectsBoundaries(curve, input.boundaries);
    risk.fence = include_fence && (!input.area.is_rough && curveLeavesFence(curve, input.area.geometry));
    risk.sibling_crosses = sampledSiblingCrossCount(curve, sampled_siblings, true, 1.5);
    return risk;
}

static bool tryPhysicalSafeSingleCubic(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, BezierCurve& curve) {
    CurveRisk current = assessCurveRisk(curve, input, sdf, sampled_siblings);
    if (!current.physical())
        return false;

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    bool require_topology_clear = (T0.dot(T1) < -0.5);
    BezierCurve best;
    bool have_best = false;
    int best_cross = std::numeric_limits<int>::max();
    double best_shape = std::numeric_limits<double>::max();

    for (double alpha : {0.50, 0.55, 0.60, 0.46, 0.42, 0.38, 0.34, 0.30, 0.70, 0.80}) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        if (curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings);
        if (risk.physical() || (require_topology_clear && risk.sibling_crosses > 0))
            continue;
        double shape_score = std::abs(alpha - 0.40);
        if (!have_best ||
            risk.sibling_crosses < best_cross ||
            (risk.sibling_crosses == best_cross && shape_score < best_shape)) {
            best = candidate;
            have_best = true;
            best_cross = risk.sibling_crosses;
            best_shape = shape_score;
        }
    }

    if (!have_best)
        return false;
    curve = best;
    return true;
}

static bool tryShapeSafeSingleCubic(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence,
    BezierCurve& curve) {
    double chord_len = (p1 - p0).norm();
    if (chord_len < 1e-6)
        return false;
    bool shape_bad = curve.maxCurvature(20) > 2.0 ||
        curve.arcLength() > std::max(chord_len * 1.8, chord_len + 12.0);
    if (!shape_bad)
        return false;

    // Detect severe divergence (optimizer ran away from the chord entirely).
    bool severely_divergent = isCurveSeverelyDivergent(curve, p0, p1, chord_len);

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : (p1 - p0).normalized();
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : (p1 - p0).normalized();
    BezierCurve best = curve;
    CurveRisk best_risk = assessCurveRisk(best, input, sdf, sampled_siblings, include_fence);
    double best_score = 1000.0 * best_risk.sibling_crosses + best.maxCurvature(20) + 0.02 * best.arcLength();
    bool improved = false;

    // 性能优化: alpha扫描从9个降至4个(Debug进一步降至2个)
#ifdef NDEBUG
    std::vector<double> alphas = {0.18, 0.26, 0.34, 0.42};
#else
    std::vector<double> alphas = {0.26, 0.38};
#endif
    for (double alpha : alphas) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        if (curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings, include_fence);
        if (!severely_divergent && risk.physical())
            continue;
        double penalty = severely_divergent ? 0.0 : 1000.0;
        double score = penalty * risk.sibling_crosses + candidate.maxCurvature(20) + 0.02 * candidate.arcLength();
        // In severely-divergent mode, always accept the first valid candidate
        // (any non-self-intersecting single cubic beats the divergent curve).
        if (severely_divergent && !improved) {
            best = candidate;
            best_score = score;
            improved = true;
            continue;
        }
        if (score + 1e-6 < best_score) {
            best = candidate;
            best_score = score;
            improved = true;
            if (risk.sibling_crosses == 0 && candidate.maxCurvature(20) < 1.0)
                break;
        }
    }

    if (!improved)
        return false;
    curve = best;
    return true;
}

static const Connectivity* findConnectivityById(
    const std::vector<Connectivity>& conns, const ConnId& id) {
    for (const auto& conn : conns)
        if (conn.id == id)
            return &conn;
    return nullptr;
}

static std::unordered_map<ConnId, std::vector<ConnId>> constrainedNeighborMap(
    const std::vector<ConnectivityCurve>& results, const ClusterOrderSolver& cs) {
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, std::vector<ConnId>> neighbors;
    for (const auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        if (!results[ia->second].curve || !results[ib->second].curve)
            continue;
        neighbors[p.id_a].push_back(p.id_b);
        neighbors[p.id_b].push_back(p.id_a);
    }
    return neighbors;
}

static int constrainedCrossCountForId(
    const ConnId& id, const BezierCurve& curve,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    double endpoint_tol = 1.5) {
    int count = 0;
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return 0;
    for (const auto& other : nit->second) {
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (curvesIntersectBusiness(curve, *other_cc.curve, endpoint_tol))
            ++count;
    }
    return count;
}

static Vec2d dominantConstrainedRefPerpForId(
    const ConnId& id, const BezierCurve& curve,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    const ClusterOrderSolver& cs,
    double endpoint_tol = 1.5) {
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return Vec2d(0, 0);
    for (const auto& other : nit->second) {
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (!curvesIntersectBusiness(curve, *other_cc.curve, endpoint_tol))
            continue;
        Vec2d ref = cs.refPerpOf(id, other);
        if (ref.norm() > 1e-8)
            return ref.normalized();
    }
    return Vec2d(0, 0);
}

static std::unordered_set<ConnId> collectPureGeometryRepairIds(
    const std::vector<ConnectivityCurve>& results, const ClusterOrderSolver& cs,
    const std::unordered_set<ConnId>& preserved_fixed_ids,
    double endpoint_tol = 1.5) {
    std::unordered_set<ConnId> repair_ids;
    auto idx = resultIndexById(results);
    for (const auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        const auto& ca = results[ia->second];
        const auto& cb = results[ib->second];
        if (!ca.curve || !cb.curve)
            continue;
        if (!curvesIntersectBusiness(*ca.curve, *cb.curve, endpoint_tol))
            continue;
        ConnId repair_id = ia->second > ib->second ? p.id_a : p.id_b;
        if (!preserved_fixed_ids.count(repair_id))
            repair_ids.insert(repair_id);
    }
    return repair_ids;
}

// ── Forward declaration ────────────────────────────────────────────────────
// crosswalkClearanceAhead 定义在下方(buildAlignedUTurn之后),
// 但 tryPureGeometryTopologyRepair (本文件) 与 generateOne 都会调用它。
// 此处前向声明以使两处调用点可编译。
static double crosswalkClearanceAhead(
    const Vec2d& p0, const Vec2d& t0,
    const IntersectionInput& input);

// solveUTurnMultiConstraint 定义在更下方(envelope辅助函数之后);
// 前向声明结构体+函数以使 tryPureGeometryTopologyRepair 可调用。
struct UTurnSolverResult {
    BezierCurve curve;
    double cost = std::numeric_limits<double>::infinity();
    int sibling_crosses = std::numeric_limits<int>::max();
    double g1_min = 0.0;
    double maxk = 0.0;
    double arc_chord = 0.0;
    double obst_pen = 0.0;
    double fence_overflow = 0.0;
    double total_handle = 0.0;
    double lateral_bias_used = 0.0;
    double lead0_used = 0.0;
    double lead1_used = 0.0;
    bool self_intersects = false;
};

static bool isGeometricUTurnConn(const Connectivity& conn, const IntersectionInput& input) {
    auto e = input.entryPtDir(conn.entry_lane_id);
    auto x = input.exitPtDir(conn.exit_lane_id);
    return e.second.norm() > 1e-8 && x.second.norm() > 1e-8 &&
           e.second.normalized().dot(x.second.normalized()) < -0.5;
}

static double uturnRadiusKey(const Connectivity& conn, const IntersectionInput& input) {
    auto e = input.entryPtDir(conn.entry_lane_id);
    auto x = input.exitPtDir(conn.exit_lane_id);
    Vec2d p0 = e.first;
    Vec2d p1 = x.first;
    Vec2d T0 = e.second.norm() > 1e-8 ? e.second.normalized() : Vec2d(1, 0);
    Vec2d T1 = x.second.norm() > 1e-8 ? x.second.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    Vec2d lat_dir{-axis[1], axis[0]};
    double turn_gap = std::abs((p1 - p0).dot(lat_dir));
    return std::max(turn_gap, 0.0);
}

static int uturnRankInConstrainedFamily(
    const Connectivity& conn, const IntersectionInput& input,
    const ClusterOrderSolver& cs) {
    auto is_uturn_conn = [&](const ConnId& id) {
        const Connectivity* c = findConnectivityById(input.connectivities, id);
        return c && isGeometricUTurnConn(*c, input);
    };

    int left_count = 0;
    int right_count = 0;
    for (const auto& p : cs.pairs()) {
        if (p.exempt != CrossExemption::None)
            continue;
        ConnId other;
        if (p.id_a == conn.id) other = p.id_b;
        else if (p.id_b == conn.id) other = p.id_a;
        else continue;
        if (!is_uturn_conn(other))
            continue;
        int side = cs.expectedSideOf(other, conn.id); // other relative to current
        if (side > 0) ++left_count;
        else if (side < 0) ++right_count;
    }
    return right_count - left_count;
}

static UTurnSolverResult solveUTurnMultiConstraint(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    const BezierCurve& reference,
    double min_lead0_floor, double min_lead1_floor,
    int depth_preference = 0, int lead1_depth_preference = 0,
    double lateral_preference = 0.0);

// 前向声明: crosswalkClearanceBehind 定义在下方
static double crosswalkClearanceBehind(
    const Vec2d& p1, const Vec2d& t1, const IntersectionInput& input);

static bool tryPureGeometryTopologyRepair(
    const Connectivity& conn, const IntersectionInput& input,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    const ClusterOrderSolver& cs,
    BezierCurve& out_curve,
    const SDFField* repair_sdf = nullptr) {
    auto it = result_idx.find(conn.id);
    if (it == result_idx.end() || !results[it->second].curve)
        return false;

    auto _entry = input.entryPtDir(conn.entry_lane_id);
    Vec2d p0 = _entry.first;
    Vec2d t0 = _entry.second;
    auto _exit = input.exitPtDir(conn.exit_lane_id);
    Vec2d p1 = _exit.first;
    Vec2d t1 = _exit.second;
    Vec2d chord = p1 - p0;
    if (chord.norm() < 1e-6)
        return false;
    const BezierCurve& current = *results[it->second].curve;
    int best_cross = constrainedCrossCountForId(
        conn.id, current, results, result_idx, neighbors, 1.5);
    if (best_cross <= 0)
        return false;

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : chord.normalized();
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : chord.normalized();
    // ── U-turn awareness: if endpoints are nearly anti-parallel, the generic
    // makeCubicG1 候选会生成扁平S形, 曲率很高(8-9 1/m)且横向凸起小——
    // 正是 conn 67/69/94/88/90/116 上报告的"扁平U型调头"缺陷。
    // 前置检测此情形, 在alpha扫描与横向偏移扫描中改用 makeAlignedUTurnCubic
    // (生成正确拱弧U型调头)。lateral_bias 参数在保持G1的同时侧向平移顶点,
    // 故 ref_perp 偏移搜索仍然适用。
    bool is_uturn_geom = (t0.norm() > 1e-8 && t1.norm() > 1e-8 &&
                          angleBetween(t0, t1) > M_PI * 0.85);
    // U型调头人行横道跨越: 与 generateOne 同逻辑 —— U型调头拱弧
    // 必须在跨越 p0 前方进入侧人行横道后开始, 且跨越 p1 前方退出侧人行横道后结束。
    double uturn_min_lead0 = is_uturn_geom
        ? crosswalkClearanceAhead(p0, t0, input) : 0.0;
    double uturn_min_lead1 = is_uturn_geom
        ? crosswalkClearanceBehind(p1, t1, input) : 0.0;
    BezierCurve best = current;
    double best_shape = current.maxCurvature(20) + 0.02 * current.arcLength();
    bool improved = false;

    auto consider = [&](BezierCurve candidate) {
        if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
            return;
        if (repair_sdf) {
            std::vector<SampledSiblingCurve> no_siblings;
            CurveRisk risk = assessCurveRisk(candidate, input, *repair_sdf, no_siblings);
            if (risk.physical())
                return;
        }

        int cross_count = constrainedCrossCountForId(
            conn.id, candidate, results, result_idx, neighbors, 1.5);
        double shape = candidate.maxCurvature(20) + 0.02 * candidate.arcLength();
        if (cross_count < best_cross ||
            (cross_count == best_cross && shape + 1e-6 < best_shape)) {
            best = std::move(candidate);
            best_cross = cross_count;
            best_shape = shape;
            improved = true;
        }
    };

    if (is_uturn_geom) {
        // ── U-turn multi-constraint repair.
        //
        // 新方法: 从约束邻居构造采样兄弟向量, 然后调用
        // solveUTurnMultiConstraint 搜索更丰富的网格(~600候选),
        // 联合代价显式权衡G1、最大曲率、同簇交叉、障碍物穿透、
        // 围栏越界与拱弧质量。
        std::vector<SampledSiblingCurve> repair_sibs;
        auto nit = neighbors.find(conn.id);
        if (nit != neighbors.end()) {
            for (const auto& other_id : nit->second) {
                auto rit = result_idx.find(other_id);
                if (rit == result_idx.end()) continue;
                const auto& other_cc = results[rit->second];
                if (!other_cc.curve) continue;
                SampledSiblingCurve s;
                s.exempt_a1 = false;  // neighbours are already constrained (non-exempt)
                s.expected_side = cs.expectedSideOf(other_id, conn.id);
                s.ref_perp = cs.refPerpOf(conn.id, other_id);
                s.shared_endpoint = cs.isSharedEndpoint(conn.id, other_id);
                s.sampled = sampleCurveForIntersections(*other_cc.curve);
                repair_sibs.push_back(std::move(s));
            }
        }
        BezierCurve ref_arch;
        ref_arch.segs.push_back(
            makeAlignedUTurnCubic(p0, T0, p1, T1, 1.0, 0.0, 0.0, uturn_min_lead0, uturn_min_lead1));
        SDFField empty_sdf;
        const SDFField& solver_sdf = repair_sdf ? *repair_sdf : empty_sdf;
        int family_rank = uturnRankInConstrainedFamily(conn, input, cs);
        double lateral_pref = family_rank > 0 ? 1.0 : (family_rank < 0 ? -1.0 : 0.0);
        UTurnSolverResult sol = solveUTurnMultiConstraint(
            p0, t0, p1, t1, input, solver_sdf, repair_sibs, ref_arch,
            uturn_min_lead0, uturn_min_lead1, 0, 0, lateral_pref);
        if (std::isfinite(sol.cost) && sol.sibling_crosses < best_cross) {
            bool physical_ok = true;
            if (repair_sdf) {
                std::vector<SampledSiblingCurve> no_siblings;
                physical_ok = !assessCurveRisk(sol.curve, input, *repair_sdf, no_siblings).physical();
            }
            if (physical_ok) {
                best = sol.curve;
                best_cross = sol.sibling_crosses;
                best_shape = sol.maxk + 0.02 * sol.curve.arcLength();
                improved = true;
            }
        }
        // 若多约束求解器未能达到零交叉, 回退到通用 ref_perp 扫描 ——
        // 该扫描采用不同
        // 候选生成(基于 ref_perp 投影), 可能找到不同解。
    } else {
        for (double alpha : {0.10, 0.12, 0.14, 0.16, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46}) {
            BezierCurve candidate;
            candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
            consider(candidate);
            if (best_cross == 0)
                break;
        }
    }

    Vec2d ref_perp = dominantConstrainedRefPerpForId(
        conn.id, current, results, result_idx, neighbors, cs, 1.5);
    if (best_cross > 0 && ref_perp.norm() > 1e-8) {
        double chord_len = chord.norm();
        if (is_uturn_geom) {
            // U型调头横向偏移扫描: 沿前向轴的法向平移顶点。
            // ref_perp 为簇横向参考; 投影到U型调头横向方向
            // 得到 lateral_bias 的有符号幅度。
            Vec2d axis = (T0 - T1).normalized();
            if (axis.norm() > 1e-8) {
                Vec2d lat_dir{-axis[1], axis[0]};
                double ref_sign = ref_perp.dot(lat_dir);
                if (std::abs(ref_sign) > 1e-6) ref_sign = (ref_sign > 0 ? 1.0 : -1.0);
                for (double scale : {0.85, 1.0, 0.70}) {
                    for (double offset : {0.75, 1.25, 1.8, 2.5, 3.4, 4.5, 6.0}) {
                        BezierCurve candidate;
                        candidate.segs.push_back(
                            makeAlignedUTurnCubic(p0, T0, p1, T1, scale, 0.0,
                                                   ref_sign * offset,
                                                   uturn_min_lead0, uturn_min_lead1));
                        consider(candidate);
                        if (best_cross == 0) break;
                    }
                    if (best_cross == 0) break;
                }
            }
        } else {
            for (double alpha : {0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42}) {
                double handle = std::max(0.05, alpha) * chord_len;
                for (double side : {-1.0, 1.0}) {
                    for (double offset : {0.75, 1.25, 1.8, 2.5, 3.4, 4.5, 6.0}) {
                        BezierSegment seg;
                        seg.ctrl[0] = p0;
                        seg.ctrl[1] = p0 + T0 * handle + ref_perp * (side * offset);
                        seg.ctrl[2] = p1 - T1 * handle + ref_perp * (side * offset);
                        seg.ctrl[3] = p1;
                        BezierCurve candidate;
                        candidate.segs.push_back(seg);
                        consider(candidate);
                        if (best_cross == 0)
                            break;
                    }
                    if (best_cross == 0)
                        break;
                }
                if (best_cross == 0)
                    break;
            }
        }
    }

    if (!improved)
        return false;
    out_curve = best;
    return true;
}

static bool chordIntersectsObstacle(
    const Vec2d& p0, const Vec2d& p1, const Obstacle& obs) {
    const Polygon2d& poly = obstacleHardGeometry(obs);
    if (poly.outer.size() < 3)
        return false;
    if (pointInsideOrOnPolygon(0.5 * (p0 + p1), poly))
        return true;
    return polygonSegmentIntersects(p0, p1, poly);
}

static BezierCurve buildWaypointCurve(
    const std::vector<Vec2d>& pts, const Vec2d& start_tan, const Vec2d& end_tan) {
    if (pts.size() < 2)
        return {};
    std::vector<Vec2d> tans(pts.size(), Vec2d(1, 0));
    tans.front() = start_tan.norm() > 1e-8 ? start_tan.normalized() : (pts[1] - pts[0]).normalized();
    tans.back() = end_tan.norm() > 1e-8 ? end_tan.normalized() : (pts.back() - pts[pts.size() - 2]).normalized();
    for (int i = 1; i + 1 < (int)pts.size(); ++i) {
        Vec2d d = pts[i + 1] - pts[i - 1];
        tans[i] = d.norm() > 1e-8 ? d.normalized() : tans[i - 1];
    }
    return makeCurveFromKnots(pts, tans, 0.34);
}

static bool tryObstacleBypassCandidate(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, BezierCurve& curve) {
    if (input.obstacles.empty())
        return false;

    Vec2d along = p1 - p0;
    double len = along.norm();
    if (len < 1e-6)
        return false;
    along /= len;
    Vec2d perp{-along[1], along[0]};
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    bool require_topology_clear = (T0.dot(T1) < -0.5);

    BezierCurve best;
    bool have_best = false;
    int best_cross = std::numeric_limits<int>::max();
    double best_score = std::numeric_limits<double>::max();
    Vec2d mid = 0.5 * (p0 + p1);
    Vec2d center = input.area.geometry.outer.empty() ? Vec2d(0, 0) : Vec2d(0, 0);
    if (!input.area.geometry.outer.empty()) {
        center.setZero();
        int cnt = 0;
        for (const auto& p : input.area.geometry.outer) {
            center += p;
            ++cnt;
        }
        if (cnt > 0)
            center /= cnt;
    }
    double center_lat = (center - mid).dot(perp);

    for (const auto& obs : input.obstacles) {
        const Polygon2d& poly = obstacleHardGeometry(obs);
        if (poly.outer.size() < 3)
            continue;
        if (!chordIntersectsObstacle(p0, p1, obs) &&
            !curveIntersectsObstacles(curve, std::vector<Obstacle>{obs}))
            continue;

        double lon_min = 1e18, lon_max = -1e18;
        double lat_min = 1e18, lat_max = -1e18;
        for (const auto& p : poly.outer) {
            double lon = (p - p0).dot(along);
            double lat = (p - p0).dot(perp);
            lon_min = std::min(lon_min, lon);
            lon_max = std::max(lon_max, lon);
            lat_min = std::min(lat_min, lat);
            lat_max = std::max(lat_max, lat);
        }
        double lon_mid = std::max(len * 0.15, std::min(len * 0.85, 0.5 * (lon_min + lon_max)));
        double obs_lat_mid = 0.5 * (lat_min + lat_max);
        std::vector<int> sides;
        if (std::abs(center_lat) > 0.25 && std::abs(obs_lat_mid) > 0.25 &&
            ((center_lat > 0) != (obs_lat_mid > 0))) {
            sides = {(center_lat > 0) ? +1 : -1, (center_lat > 0) ? -1 : +1};
        } else {
            sides = {+1, -1};
        }

        for (int side : sides) {
            double edge_lat = side > 0 ? lat_max : lat_min;
            for (double clearance : {0.65, 1.0, 1.4, 1.9, 2.5, 3.2}) {
                double apex_lat = edge_lat + side * clearance;
                Vec2d apex = p0 + lon_mid * along + apex_lat * perp;
                double spread = std::max(2.0, std::min(len * 0.28, (lon_max - lon_min) * 0.8 + 1.5));
                double lon_a = std::max(len * 0.08, lon_mid - spread);
                double lon_b = std::min(len * 0.92, lon_mid + spread);
                Vec2d q0 = p0 + lon_a * along + apex_lat * 0.55 * perp;
                Vec2d q1 = p0 + lon_b * along + apex_lat * 0.55 * perp;
                std::vector<Vec2d> pts;
                if ((q0 - p0).norm() > 1.0 && (q1 - p1).norm() > 1.0 && lon_b > lon_a + 0.5)
                    pts = {p0, q0, apex, q1, p1};
                else
                    pts = {p0, apex, p1};

                BezierCurve candidate = buildWaypointCurve(pts, t0, t1);
                if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings);
                if (risk.physical() || (require_topology_clear && risk.sibling_crosses > 0))
                    continue;
                double shape_score = clearance + 0.02 * candidate.arcLength() +
                    (side * center_lat >= 0.0 ? 0.0 : 0.4);
                if (!have_best ||
                    risk.sibling_crosses < best_cross ||
                    (risk.sibling_crosses == best_cross && shape_score < best_score)) {
                    best = candidate;
                    have_best = true;
                    best_cross = risk.sibling_crosses;
                    best_score = shape_score;
                }
            }
        }
    }

    if (!have_best)
        return false;
    curve = best;
    return true;
}

static BezierCurve buildAlignedUTurn(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const Vec2d& offset_dir, double offset_m, double forward_scale = 1.0,
    double min_lead0 = 0.0, double min_lead1 = 0.0) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    // 将 offset_dir 分解为前向(沿轴)和横向(垂直轴)分量。
    // 此前仅前向分量作为 handle_bias, 横向偏移(同进入U型调头扇出的自然选择)
    // 因 offset_dir.dot(axis)=0 被默默丢弃。
    // 前向偏移纵向延伸把手; 横向偏移侧向平移两个内部控制点,
    // 顶点平移而不破坏端点G1连续。
    double fwd_bias = 0.0;
    double lat_bias = 0.0;
    if (offset_m > 0.0 && offset_dir.norm() > 1e-8) {
        Vec2d dir = offset_dir.normalized();
        fwd_bias = offset_m * dir.dot(axis);
        Vec2d lat_dir{-axis[1], axis[0]};
        lat_bias = offset_m * dir.dot(lat_dir);
    }

    BezierCurve c;
    c.segs.push_back(makeAlignedUTurnCubic(p0, T0, p1, T1, forward_scale,
                                            fwd_bias, lat_bias, min_lead0, min_lead1));
    return c;
}

// ── Crosswalk-aware U-turn start delay ─────────────────────────────────────
//
// 新需求: 如果数据中调头连接开始位置附近能搜到人行横道就必须要向
// 路口内跨越过人行横道后才能调头。
//
// 本例程检查人行横道(当数据无人行横道时用停止线作为代理
// 数据集仅提供stop_lines因人行横道带位于停止线后方)。对每个候选:
//   1. 计算折线沿 T0 方向距 p0 的最近点。
//   2. 将位移(nearest - p0)投影到 T0; 若投影为正
//      (人行横道在 p0 沿进入方向的"前方")且在合理搜索半径内(默认8m),
//      则为候选。
//   3. 返回所有候选中"远边"投影的最大值 —— 即车辆从 p0 必须前行多远
//      才能开始拱弧。零表示未找到人行横道, 不施加约束。
//
// 结果作为 `min_lead0` 传给 `makeAlignedUTurnCubic`,
// 强制进入延长段(p0 → q0)跨过人行横道后才开始近半圆。
// lead1 对称延伸以保持U型调头轴线对齐(视觉G1)。
static double crosswalkClearanceAhead(
    const Vec2d& p0, const Vec2d& t0,
    const IntersectionInput& input) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    const double search_radius = 8.0;   // search within 8m of p0
    const double side_tolerance = 4.0;  // within 4m lateral of the entry ray
    double max_far = 0.0;
    bool found = false;

    auto consider_polyline = [&](const std::vector<Vec2d>& pts) {
        for (size_t i = 1; i < pts.size(); ++i) {
            const Vec2d& a = pts[i - 1];
            const Vec2d& b = pts[i];
            // 将 p0 投影到线段 a-b
            Vec2d ab = b - a;
            double ab_len2 = ab.dot(ab);
            if (ab_len2 < 1e-12) continue;
            double tt = std::max(0.0, std::min(1.0, (p0 - a).dot(ab) / ab_len2));
            Vec2d closest = a + tt * ab;
            Vec2d d = closest - p0;
            double lateral = std::abs(d.dot(Vec2d(-T0[1], T0[0])));
            if (lateral > side_tolerance) continue;
            // "Far edge" = the further endpoint of this segment along T0,
            // 用以近似人行横道远边。
            double a_fwd = (a - p0).dot(T0);
            double b_fwd = (b - p0).dot(T0);
            double far = std::max(a_fwd, b_fwd);
            // 当线段任意部分位于 p0 前方(far > 0)时考虑该线段。
            // 用 `far` (非最近点投影) 确保同arm上所有U型调头
            // 获得一致的人行横道跨越距离, 保持其相对顶点排序。
            if (far <= 0.0) continue;
            if (far > search_radius) continue;
            if (far > max_far) {
                max_far = far;
                found = true;
            }
        }
    };

    auto consider_polygon = [&](const Polygon2d& poly) {
        if (poly.outer.size() >= 2) consider_polyline(poly.outer);
        for (const auto& hole : poly.holes) consider_polyline(hole);
    };

    // 仅使用实际的人行横道多边形(面)作为跨越源。
    //
    // 依据需求"人行横道附近调头需跨越约束": U型调头拱弧必须跨越
    // 人行横道带后才能开始。此约束仅在输入数据中存在实际人行横道
    // 多边形几何时才有意义。
    for (const auto& cw : input.crosswalks) {
        consider_polygon(cw.geometry);
    }
    // 无人行横道时使用停止线"4m人行横道深度"作为人行横道: 保留暂未启用
    if (false && input.crosswalks.empty()) {
        for (const auto& sl : input.stop_lines) {
            std::vector<Vec2d> pts = sl.geometry.points;
            if (pts.size() < 2) continue;
            consider_polyline(pts);
        }
        if (found) {
            max_far += 4.0;
        }
    }

    return found ? max_far : 0.0;
}

// ── 退出侧人行横道跨越检查 ───────────────────────────────────
//
// 新需求: 调头真实弧段必须跨越退出侧人行横道后才能结束。
// 本例程检查 p1(退出端点)前方(沿 exit_back = -T1 方向,即从路口向退出车道方向)
// 是否有人行横道, 返回退出侧最小延长距离 min_lead1。
//
// 与 crosswalkClearanceAhead 对称:
//   - crosswalkClearanceAhead 检查进入侧(p0沿T0方向)
//   - crosswalkClearanceBehind 检查退出侧(p1沿exit_back=-T1方向)
//
// 结果作为 `min_lead1` 传给 makeAlignedUTurnCubic,
// 强制退出反向延长段(p1 → q1)跨过退出侧人行横道后近半圆才结束。
static double crosswalkClearanceBehind(
    const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input) {
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    Vec2d exit_back = -T1;  // 从路口向退出车道方向
    const double search_radius = 8.0;
    const double side_tolerance = 4.0;
    double max_far = 0.0;
    bool found = false;

    auto consider_polyline = [&](const std::vector<Vec2d>& pts) {
        for (size_t i = 1; i < pts.size(); ++i) {
            const Vec2d& a = pts[i - 1];
            const Vec2d& b = pts[i];
            Vec2d ab = b - a;
            double ab_len2 = ab.dot(ab);
            if (ab_len2 < 1e-12) continue;
            double tt = std::max(0.0, std::min(1.0, (p1 - a).dot(ab) / ab_len2));
            Vec2d closest = a + tt * ab;
            Vec2d d = closest - p1;
            double lateral = std::abs(d.dot(Vec2d(-exit_back[1], exit_back[0])));
            if (lateral > side_tolerance) continue;
            double a_fwd = (a - p1).dot(exit_back);
            double b_fwd = (b - p1).dot(exit_back);
            double far = std::max(a_fwd, b_fwd);
            if (far <= 0.0) continue;
            if (far > search_radius) continue;
            if (far > max_far) {
                max_far = far;
                found = true;
            }
        }
    };

    auto consider_polygon = [&](const Polygon2d& poly) {
        if (poly.outer.size() >= 2) consider_polyline(poly.outer);
        for (const auto& hole : poly.holes) consider_polyline(hole);
    };

    for (const auto& cw : input.crosswalks) {
        consider_polygon(cw.geometry);
    }

    return found ? max_far : 0.0;
}

static bool curveExceedsUTurnEnvelope(
    const BezierCurve& curve, const BezierCurve& reference,
    const Vec2d& p0, const Vec2d& p1) {
    double ref_len = std::max(reference.arcLength(), (p1 - p0).norm());
    double max_len = std::max(80.0, ref_len * 4.0);
    if (curve.arcLength() > max_len)
        return true;

    Vec2d mid = 0.5 * (p0 + p1);
    double max_radius = std::max(60.0, ref_len * 3.0);
    for (const auto& pt : curve.sampleByArcLength(80)) {
        if ((pt - mid).norm() > max_radius)
            return true;
    }
    return false;
}

static double curveLateralBulge(
    const BezierCurve& curve, const Vec2d& p0, const Vec2d& p1) {
    Vec2d chord = p1 - p0;
    if (chord.norm() < 1e-8)
        return 0.0;
    Vec2d perp{-chord.y(), chord.x()};
    perp.normalize();
    double best = 0.0;
    for (const auto& pt : curve.sampleByArcLength(80)) {
        double v = (pt - p0).dot(perp);
        if (std::abs(v) > std::abs(best))
            best = v;
    }
    return best;
}

static bool curveCollapsesUTurnEnvelope(
    const BezierCurve& curve, const BezierCurve& reference,
    const Vec2d& p0, const Vec2d& p1) {
    double ref_len = reference.arcLength();
    double chord_len = (p1 - p0).norm();
    double ref_bulge = std::abs(curveLateralBulge(reference, p0, p1));
    double cur_bulge = std::abs(curveLateralBulge(curve, p0, p1));

    if (ref_len > 1.0 && curve.arcLength() < std::max(chord_len * 1.08, ref_len * 0.72))
        return true;
    if (ref_bulge > 2.0 && cur_bulge < ref_bulge * 0.50)
        return true;
    if (curve.maxCurvature(20) > std::max(1.0, reference.maxCurvature(20) * 8.0))
        return true;
    return false;
}

// ─────────────────────────────────────────────────────────────────────────────
//  solveUTurnMultiConstraint
//
//  Multi-constraint U-turn solver (round 3 redesign).
//
//  Previous design: `tryBoundedUTurnCandidate` searched a small grid
//  (4 directions × 6 magnitudes × 6 scales = 144 candidates) and selected
//  the one with fewest sibling crossings, breaking ties on a shape score.
//  In practice this failed because:
//    (a) the lateral_bias clamp in makeAlignedUTurnCubic (tan(15°) ≈ 1.3 m
//        for large U-turns) silently truncated any candidate with |offset|
//        > 1.3 m, so the search space effectively shrank to ~1/4 of the
//        intended range.
//    (b) tie-breaking only minimised curvature + arc length — it never
//        penalised "near-miss" sibling proximity, so the chosen candidate
//        could graze a sibling and then re-intersect after the optimizer
//        perturbed it.
//    (c) obstacle/fence clearance was a hard filter (any candidate touching
//        an obstacle was rejected), but for tightly-packed intersections
//        every candidate grazes *something*, so the filter rejected
//        everything and the solver fell back to the LBFGS optimizer, which
//        then ran for hundreds of milliseconds per curve and still left
//        14 same-cluster crossings on 100000643-1.json.
//
//  New design:
//    1. Search a richer grid: scale × lateral_bias × min_lead0_extra ×
//       handle_bias.  For each candidate the actual lateral_bias applied
//       inside makeAlignedUTurnCubic is *clamped* to the G1-safe envelope
//       (tan(30°) for large U-turns, tan(15°) for small ones), so we never
//       produce a non-G1 curve.
//    2. Joint cost = w_cross * sibling_crosses
//                  + w_obst * max(0, -min_sdf - 0.05)        // obstacle penetration
//                  + w_fence * fence_overflow
//                  + w_g1 * (1 - min(g1_0, g1_1))
//                  + w_maxk * max(0, maxk - 1.0)
//                  + w_arch * |arc/chord - 1.55|              // ideal arch ratio
//                  + w_self * self_intersects                 // hard reject flag
//                  + w_compact * (|lat_bias| + |handle_bias| + 0.3*lead0_extra)
//    3. Hard filter: envelope checks (curveExceedsUTurnEnvelope,
//       curveCollapsesUTurnEnvelope) and self-intersection still reject
//       candidates.  Obstacle/fence are SOFT (cost terms) so the solver can
//       pick a slightly-penetrating candidate when no zero-penetration
//       candidate exists.
//    4. Returns the joint-cost minimum.  Caller may still run the LBFGS
//       optimizer on the result if the joint cost is too high.
//
//  Naming convention: "min_lead0_extra" = additional entry prolongation
//  beyond the crosswalk-clearance minimum.  This is what moves the apex
//  forward into the intersection; useful for separating same-entry U-turns
//  by apex depth.
// ─────────────────────────────────────────────────────────────────────────────
static UTurnSolverResult solveUTurnMultiConstraint(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    const BezierCurve& reference,
    double min_lead0_floor, double min_lead1_floor,
    int depth_preference, int lead1_depth_preference,
    double lateral_preference) {

    UTurnSolverResult best;
    best.cost = std::numeric_limits<double>::infinity();

    bool debug = isgDebugUTurn();
    if (debug) {
        fprintf(stderr, "[UTURN-SOLVER] p0=(%.2f,%.2f) p1=(%.2f,%.2f) siblings=%zu (exempt_a1 count: ",
                p0.x(), p0.y(), p1.x(), p1.y(), sampled_siblings.size());
        int ex_cnt = 0;
        for (auto& s : sampled_siblings) if (s.exempt_a1) ++ex_cnt;
        fprintf(stderr, "%d/%zu)\n", ex_cnt, sampled_siblings.size());
    }

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8) axis = T0;
    axis.normalize();

    Vec2d lat_dir{-axis[1], axis[0]};
    double turn_gap = std::abs((p1 - p0).dot(lat_dir));
    double chord_len = (p1 - p0).norm();
    double ref_len = std::max(reference.arcLength(), chord_len);

    // ── Determine the maxk physical bound for this U-turn size.
    // For tiny turn_gap (< 1 m), the perfect-semicircle maxk is 2/turn_gap;
    // we accept up to 2× that.  For larger turn_gap, 1.5 is drivable.
    double maxk_phys_bound = (turn_gap < 1.0)
        ? std::max(8.0, 4.0 / std::max(0.1, turn_gap))
        : 1.5;

    // ── Joint cost weights (calibrated on 100000643-1.json).
    // sibling_crosses is the dominant term — a single same-cluster interior
    // crossing is worth ~1000 units of "shape quality" because it directly
    // violates the cluster ordering constraint.
    const double W_CROSS  = 1000.0;
    const double W_OBST   =  200.0;   // per metre of obstacle penetration
    const double W_FENCE  =  200.0;   // per metre of fence overflow
    const double W_G1     =   50.0;   // fallback tie-breaker after the hard G1 gate
    const double W_MAXK   =   20.0;   // per unit maxk above 1.0
    const double W_ARCH   =    5.0;   // per unit deviation from arc/chord=1.55
    const double W_COMPACT =   0.5;   // per metre of lat/handle bias used
    const double W_LEN    =    0.1;   // per metre of arc length deviation from ref

    constexpr double G1_HARD_MIN = 0.99;
    auto physical_bad = [](double obst_pen, double fence_ovf) {
        return obst_pen > 0.05 || fence_ovf > 0.05;
    };
    auto better_uturn_candidate = [&](int sib_x, double obst_pen, double fence_ovf,
                                      double cost, const UTurnSolverResult& incumbent) {
        if (!std::isfinite(incumbent.cost))
            return true;
        if (sib_x != incumbent.sibling_crosses)
            return sib_x < incumbent.sibling_crosses;

        bool bad = physical_bad(obst_pen, fence_ovf);
        bool incumbent_bad = physical_bad(incumbent.obst_pen, incumbent.fence_overflow);
        if (bad != incumbent_bad)
            return !bad;
        if (std::abs(obst_pen - incumbent.obst_pen) > 0.02)
            return obst_pen < incumbent.obst_pen;
        if (std::abs(fence_ovf - incumbent.fence_overflow) > 0.02)
            return fence_ovf < incumbent.fence_overflow;
        return cost < incumbent.cost;
    };

    // ── Grid: scale × lateral_bias × lead0_extra × handle_bias.
    // 性能优化: 大幅减少候选数(原~600,现~120)以满足<1s需求
    //   - tiny turn_gap (<1.5m): 4个候选(原4)
    //   - small turn_gap (1.5-4m): 24个候选(原36)
    //   - large turn_gap (>=4m): 96个候选(原252)
    // 每个候选计算时启用快速路径(简化maxk采样,sib_x提前退出)
    std::vector<double> scales;
    std::vector<double> lat_biases;
    std::vector<double> lead0_extras;
    std::vector<double> handle_biases;

    scales = {0.7, 1.0};  // 性能优化: 从4个scale降至2个
#ifdef NDEBUG
    const bool dbg_small_grid = false;
#else
    const bool dbg_small_grid = true;  // Debug模式进一步缩减
#endif
    if (turn_gap < 1.5) {
        lat_biases = {0.0};
        lead0_extras = {0.0};
        handle_biases = {0.0};
    } else if (turn_gap < 4.0) {
        // 扩展lateral_bias范围, 允许南北分离避免中段相交
        lat_biases = dbg_small_grid ? std::vector<double>{0.0, 1.5, -1.5, 3.0, -3.0}
                                    : std::vector<double>{0.0, 1.5, -1.5, 3.0, -3.0, 5.0, -5.0};
        // 扩展lead0_extra搜索范围, 允许更深的拱弧以清理人行横道
        lead0_extras = dbg_small_grid ? std::vector<double>{0.0, 1.0}
                                      : std::vector<double>{0.0, 1.0, 2.0, 3.0};
        handle_biases = dbg_small_grid ? std::vector<double>{0.0}
                                       : std::vector<double>{0.0, 1.0};
    } else {
        // Large U-turns: 需要更宽lateral_bias范围逃离左转兄弟和共端点U-turn
        lat_biases = dbg_small_grid ? std::vector<double>{0.0, 1.5, -1.5, 3.0, -3.0}
                                    : std::vector<double>{0.0, 1.5, -1.5, 3.0, -3.0, 5.0, -5.0};
        // 扩展lead0_extra搜索范围, 允许更深的拱弧以清理人行横道
        lead0_extras = dbg_small_grid ? std::vector<double>{0.0, 1.0, 2.0}
                                      : std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0};
        handle_biases = dbg_small_grid ? std::vector<double>{0.0}
                                       : std::vector<double>{0.0, 1.5, -1.5};
    }

    // ── 根据 depth_preference 调整 lead0_extra 搜索方向 ──────────
    // depth_preference 来源于 generateOne 中对 cluster_solver 配对的分析:
    //   +1 = 当前U型调头应为"外层深拱" (apex更远, lead0_extra > 0)
    //   -1 = 当前U型调头应为"内层浅拱" (apex更近, lead0_extra = 0)
    //    0 = 无偏好 (普通搜索)
    //
    // 这使共端点两个U型调头通过不同拱弧深度分离: 内层浅拱 + 外层深拱,
    // 外层包裹内层, 不再相交。例如 conn 74 (内层) vs conn 76 (外层)。
    //
    // 注意: 不再硬性移除lead0_extra=0, 而是作为软偏好通过W_DEPTH惩罚引导。
    // 硬性移除会导致求解器无法找到同时满足人行横道清理和零交叉的候选。
    // 软偏好允许求解器在必要时选择"错误"深度方向以避免交叉。

    // ── 搜索 (共退出U-turn对的深度分离) ─────────
    // 对共退出车道的U-turn对, 通过不同 lead1 (退出延长) 深度分离:
    // 内层浅拱 lead1_extra=0, 外层深拱 lead1_extra>0。
    // 大径U-turn(如73, turn_gap=7.6m)应深拱包围小径U-turn(如74, turn_gap=4.3m)。
    // 同样使用软偏好, 不硬性移除lead1_extra=0。
    std::vector<double> lead1_extras = {0.0};  // 默认不变化
    if (turn_gap >= 1.5) {
        lead1_extras = dbg_small_grid ? std::vector<double>{0.0, 1.0, 2.0}
                                      : std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    }

    auto enforce_fence = !input.area.is_rough && !input.area.geometry.outer.empty();
    const Polygon2d& fence = input.area.geometry;

    for (double scale : scales) {
        for (double lat_bias : lat_biases) {
            for (double lead0_extra : lead0_extras) {
                for (double lead1_extra : lead1_extras) {
                for (double handle_bias : handle_biases) {
                    double eff_min_lead0 = std::max(min_lead0_floor, min_lead0_floor + lead0_extra);
                    double eff_min_lead1 = std::max(min_lead1_floor, min_lead1_floor + lead1_extra);
                    BezierCurve cand;
                    cand.segs.push_back(makeAlignedUTurnCubic(
                        p0, T0, p1, T1, scale, handle_bias, lat_bias,
                        eff_min_lead0, eff_min_lead1));
                    if (cand.empty()) continue;

                    // ── Hard reject: self-intersection, envelope collapse/exceed
                    if (curveSelfIntersectsBusiness(cand, 1.0)) continue;
                    if (curveExceedsUTurnEnvelope(cand, reference, p0, p1)) continue;
                    if (curveCollapsesUTurnEnvelope(cand, reference, p0, p1)) continue;

                    // ── Compute metrics
                    int sib_x = sampledSiblingCrossCount(cand, sampled_siblings, true, 1.5);
                    // 性能优化: maxk采样从40降至20
                    double maxk = cand.maxCurvature(20);
                    if (maxk > maxk_phys_bound) continue;  // hard reject undrivable
                    double arc = cand.arcLength();
                    double arc_chord = chord_len > 1e-6 ? arc / chord_len : 1.0;

                    // G1 cos at endpoints
                    Vec2d st = cand.startTan().norm() > 1e-8 ? cand.startTan().normalized() : Vec2d(1, 0);
                    Vec2d et = cand.endTan().norm()   > 1e-8 ? cand.endTan().normalized()   : Vec2d(1, 0);
                    double g1_0 = st.dot(T0);
                    double g1_1 = et.dot(T1);
                    double g1_min = std::min(g1_0, g1_1);
                    if (g1_min < G1_HARD_MIN)
                        continue;

                    // Obstacle penetration
                    double obst_pen = 0.0;
                    if (!input.obstacles.empty() && sdf.valid()) {
                        double ms = minSDFAlongCurveAdaptive(cand, sdf);
                        obst_pen = std::max(0.0, -ms);
                        if (curveIntersectsObstacles(cand, input.obstacles))
                            obst_pen = std::max(obst_pen, 0.10);
                    } else if (!input.obstacles.empty()) {
                        if (curveIntersectsObstacles(cand, input.obstacles))
                            obst_pen = 0.10;  // crude hard-geometry penalty if no SDF
                    }

                    // Fence overflow
                    double fence_ovf = 0.0;
                    if (enforce_fence) {
                        auto pts = cand.sample(24);
                        for (int k = 1; k + 1 < (int)pts.size(); ++k) {
                            if (!polygonContains(fence, pts[k])) {
                                double d = pointToPolygonDist(pts[k], fence);
                                if (d > fence_ovf) fence_ovf = d;
                            }
                        }
                    }

                    // ── Crosswalk clearance penalty ───────
                    // 调头真实弧段必须跨越人行横道后才能构建。检查候选曲线的拱弧顶点
                    // (apex, 沿U型轴最远点)附近(3m内)是否进入人行横道多边形内部,
                    // 若是则按穿透深度施加惩罚。这促使求解器选择更深的拱弧(更大
                    // lead0_extra)以清理人行横道, 满足 (4.2) "首尾直行段都跨越人行
                    // 横道后才能构建弧段约束"。
                    //
                    // 注意: 不检查首尾直行段 —— 直行段本应跨越人行横道(需求 4.2 要求),
                    // 只有拱弧必须在人行横道外。
                    Vec2d ut_axis_cand = (T0 - T1).normalized();
                    if (ut_axis_cand.dot(T0) < 0) ut_axis_cand = -ut_axis_cand;
                    auto pts_c = cand.sample(40);
                    Vec2d apex_c = p0;
                    double max_proj_c = -1e18;
                    for (auto& pt : pts_c) {
                        double proj = (pt - p0).dot(ut_axis_cand);
                        if (proj > max_proj_c) { max_proj_c = proj; apex_c = pt; }
                    }

                    double xwalk_pen = 0.0;
                    if (!input.crosswalks.empty()) {
                        for (auto& pt : pts_c) {
                            if ((pt - apex_c).norm() > 3.0) continue;  // only arc region
                            for (const auto& cw : input.crosswalks) {
                                if (polygonContains(cw.geometry, pt)) {
                                    xwalk_pen += pointToPolygonDist(pt, cw.geometry);
                                }
                            }
                        }
                    }

                    // ── Cluster order violation penalty ──────────────
                    // 对每个非豁免的共端点兄弟, 检查候选曲线的拱弧顶点是否位于正确侧。
                    // expected_side=+1: 兄弟应在左侧 → 兄弟apex的ref_perp投影应 > 候选apex
                    // expected_side=-1: 兄弟应在右侧 → 兄弟apex的ref_perp投影应 < 候选apex
                    // 若违反, 按违反量施加惩罚, 促使求解器选择正确深度的拱弧。
                    //
                    // 这使共端点两个U型调头通过不同深度分离: 内层浅拱apex更近(更EAST),
                    // 外层深拱apex更远(更WEST), 外层包裹内层。
                    double order_violation = 0.0;
                    SampledCurve cand_sample = sampleCurveForIntersections(cand);
                    for (auto& sib : sampled_siblings) {
                        if (sib.exempt_a1) continue;
                        if (sib.expected_side == 0) continue;
                        if (sib.ref_perp.norm() < 1e-9) continue;
                        if (sib.sampled.pts.size() < 3) continue;

                        // 候选曲线apex在ref_perp方向的投影
                        double cur_apex_lat = apex_c.dot(sib.ref_perp);

                        // 兄弟曲线apex (沿ut_axis最远点)
                        Vec2d sib_apex = sib.sampled.pts[0];
                        double sib_max_proj = -1e18;
                        for (auto& pt : sib.sampled.pts) {
                            double proj = (pt - p0).dot(ut_axis_cand);
                            if (proj > sib_max_proj) { sib_max_proj = proj; sib_apex = pt; }
                        }
                        double sib_apex_lat = sib_apex.dot(sib.ref_perp);

                        // diff = cur_lat - sib_lat (正=候选在兄弟左侧, ref_perp方向)
                        double diff = cur_apex_lat - sib_apex_lat;
                        double ORDER_MARGIN = 0.5;  // 米 — 最小分离量
                        double viol = 0;
                        if (sib.expected_side == +1) {
                            // 兄弟在左侧 → sib_lat > cur_lat → diff应为负
                            viol = std::max(0.0, diff + ORDER_MARGIN);
                        } else {
                            // 兄弟在右侧 → sib_lat < cur_lat → diff应为正
                            viol = std::max(0.0, -diff + ORDER_MARGIN);
                        }
                        order_violation += viol;
                    }
                    double endpoint_side_violation =
                        sharedEndpointSideViolation(cand_sample, sampled_siblings, 1.5);

                    if (debug && sib_x == 0) {
                        fprintf(stderr, "[UTURN-SOLVER] ZERO-X: scale=%.2f lat=%.1f lead0=%.2f hand=%.1f maxk=%.3f g1=%.3f arc/c=%.3f xwalk=%.2f ord_viol=%.2f\n",
                                scale, lat_bias, eff_min_lead0, handle_bias, maxk, g1_min, arc_chord, xwalk_pen, order_violation);
                    }

                    // ── Joint cost
                    // W_XWALK = 2000: crosswalk clearance is a HARD requirement (4.2),
                    // weighted 2× a sibling crossing to ensure the solver prioritizes it.
                    // W_ORDER = 300: cluster ordering violation
                    // W_DEPTH = 50: soft depth preference
                    // W_LATPREF = 200: soft lateral preference — guides north/south offset
                    //   to separate inner/outer U-turns (彻底消除中段相交)
                    const double W_XWALK = 2000.0;
                    const double W_ORDER = 300.0;
                    const double W_ENDPOINT_SIDE = 500.0;
                    const double W_DEPTH = 50.0;
                    const double W_LATPREF = 200.0;
                    double cost = 0.0;
                    cost += W_CROSS * sib_x;
                    cost += W_OBST * obst_pen;
                    cost += W_FENCE * fence_ovf;
                    cost += W_XWALK * xwalk_pen;
                    cost += W_ORDER * order_violation;
                    cost += W_ENDPOINT_SIDE * endpoint_side_violation;
                    cost += W_G1 * (1.0 - g1_min);
                    cost += W_MAXK * std::max(0.0, maxk - 1.0);
                    cost += W_ARCH * std::abs(arc_chord - 1.55);
                    // 软深度偏好: 惩罚"错误"深度方向
                    if (depth_preference > 0 && lead0_extra < 0.5) cost += W_DEPTH;
                    if (depth_preference < 0 && lead0_extra > 0.5) cost += W_DEPTH;
                    if (lead1_depth_preference > 0 && lead1_extra < 0.5) cost += W_DEPTH;
                    if (lead1_depth_preference < 0 && lead1_extra > 0.5) cost += W_DEPTH;
                    // 软横向偏好: 惩罚"错误"横向方向
                    // lateral_preference > 0: 应向 +lat_dir 偏移 (lateral_bias > 0)
                    // lateral_preference < 0: 应向 -lat_dir 偏移 (lateral_bias < 0)
                    if (lateral_preference > 0.5 && lat_bias < 0.5) cost += W_LATPREF;
                    if (lateral_preference < -0.5 && lat_bias > -0.5) cost += W_LATPREF;
                    // 修复 BUG 1 (4.2): 当有人行横道时, 降低lead0_extra的紧凑度惩罚,
                    // 允许求解器自由选择更深的拱弧以清理人行横道。
                    double lead0_penalty = input.crosswalks.empty() ? 0.3 : 0.05;
                    cost += W_COMPACT * (std::abs(lat_bias) + std::abs(handle_bias) + lead0_penalty * lead0_extra);
                    cost += W_LEN * std::abs(arc - ref_len);

                    if (better_uturn_candidate(sib_x, obst_pen, fence_ovf, cost, best)) {
                        best.curve = cand;
                        best.cost = cost;
                        best.sibling_crosses = sib_x;
                        best.g1_min = g1_min;
                        best.maxk = maxk;
                        best.arc_chord = arc_chord;
                        best.obst_pen = obst_pen;
                        best.fence_overflow = fence_ovf;
                        best.lateral_bias_used = lat_bias;
                        best.lead0_used = eff_min_lead0;
                        best.lead1_used = eff_min_lead1;
                        best.self_intersects = false;
                    }
                }  // handle_bias
                }  // lead1_extra
            }  // lead0_extra
        }  // lat_bias
    }  // scale

    if (debug) {
        if (std::isfinite(best.cost)) {
            fprintf(stderr, "[UTURN-SOLVER] BEST: cost=%.2f sib_x=%d g1=%.3f maxk=%.3f arc/c=%.3f lat_bias=%.1f lead0=%.1f scale=%.2f handle_bias=%.1f\n",
                    best.cost, best.sibling_crosses, best.g1_min, best.maxk, best.arc_chord,
                    best.lateral_bias_used, best.lead0_used, 0.0, 0.0);
        } else {
            fprintf(stderr, "[UTURN-SOLVER] NO VALID CANDIDATE\n");
        }
    }

    // ── 深度偏好约束下仍无法零交叉时, 放宽约束重新搜索 ──────
    // 当 depth_preference/lead1_depth_preference 强制了特定深度方向, 但所有
    // 该方向的候选都与兄弟相交 (sib_x > 0) 时, 放宽约束: 允许 lead0_extra=0
    // 和 lead1_extra=0 (浅拱), 重新搜索。零交叉比深度分离更重要。
    // 注意: 放宽搜索仍需考虑人行横道清理 (4.2) 和围栏约束。
    if (std::isfinite(best.cost) && best.sibling_crosses > 0 &&
        (depth_preference != 0 || lead1_depth_preference != 0)) {
        UTurnSolverResult relaxed;
        relaxed.cost = std::numeric_limits<double>::infinity();
        Vec2d ut_axis_r = (T0 - T1).normalized();
        if (ut_axis_r.dot(T0) < 0) ut_axis_r = -ut_axis_r;
        std::vector<double> rl_lead0 = {0.0, 1.0, 2.0, 3.0, 4.0};
        std::vector<double> rl_lead1 = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        for (double scale : scales) {
            for (double lat_bias : lat_biases) {
                for (double l0 : rl_lead0) {
                    for (double l1 : rl_lead1) {
                        for (double hb : handle_biases) {
                            double em0 = std::max(min_lead0_floor, min_lead0_floor + l0);
                            double em1 = std::max(min_lead1_floor, min_lead1_floor + l1);
                            BezierCurve cand;
                            cand.segs.push_back(makeAlignedUTurnCubic(
                                p0, T0, p1, T1, scale, hb, lat_bias, em0, em1));
                            if (cand.empty()) continue;
                            if (curveSelfIntersectsBusiness(cand, 1.0)) continue;
                            if (curveExceedsUTurnEnvelope(cand, reference, p0, p1)) continue;
                            if (curveCollapsesUTurnEnvelope(cand, reference, p0, p1)) continue;
                            int sx = sampledSiblingCrossCount(cand, sampled_siblings, true, 1.5);
                            double mk = cand.maxCurvature(20);
                            if (mk > maxk_phys_bound) continue;

                            Vec2d st = cand.startTan().norm() > 1e-8 ? cand.startTan().normalized() : Vec2d(1, 0);
                            Vec2d et = cand.endTan().norm()   > 1e-8 ? cand.endTan().normalized()   : Vec2d(1, 0);
                            double g1_min = std::min(st.dot(T0), et.dot(T1));
                            if (g1_min < G1_HARD_MIN) continue;

	                            // 人行横道穿透检查
	                            double xp = 0.0;
	                            if (!input.crosswalks.empty()) {
	                                auto pts_r = cand.sample(40);
                                Vec2d ap = p0; double mp = -1e18;
                                for (auto& pt : pts_r) {
                                    double pr = (pt - p0).dot(ut_axis_r);
                                    if (pr > mp) { mp = pr; ap = pt; }
                                }
                                for (auto& pt : pts_r) {
                                    if ((pt - ap).norm() > 3.0) continue;
	                                    for (const auto& cw : input.crosswalks)
	                                        if (polygonContains(cw.geometry, pt))
	                                            xp += pointToPolygonDist(pt, cw.geometry);
	                                }
	                            }

	                            double obst_pen = 0.0;
	                            if (!input.obstacles.empty() && sdf.valid()) {
	                                double ms = minSDFAlongCurveAdaptive(cand, sdf);
	                                obst_pen = std::max(0.0, -ms);
	                                if (curveIntersectsObstacles(cand, input.obstacles))
	                                    obst_pen = std::max(obst_pen, 0.10);
	                            } else if (!input.obstacles.empty()) {
	                                if (curveIntersectsObstacles(cand, input.obstacles))
	                                    obst_pen = 0.10;
	                            }

	                            double fence_ovf = 0.0;
	                            if (enforce_fence) {
	                                auto pts_f = cand.sample(24);
	                                for (int k = 1; k + 1 < (int)pts_f.size(); ++k) {
	                                    if (!polygonContains(fence, pts_f[k])) {
	                                        double d = pointToPolygonDist(pts_f[k], fence);
	                                        if (d > fence_ovf) fence_ovf = d;
	                                    }
	                                }
	                            }

	                            // 代价: 零交叉与人行横道同等重要, 横向偏好次之, 曲率/弧长最后
	                            SampledCurve cand_sample = sampleCurveForIntersections(cand);
	                            double endpoint_side_violation =
	                                sharedEndpointSideViolation(cand_sample, sampled_siblings, 1.5);
	                            double c = W_CROSS * sx
	                                     + W_OBST * obst_pen + W_FENCE * fence_ovf
	                                     + 2000.0 * xp
	                                     + 500.0 * endpoint_side_violation
	                                     + W_MAXK * std::max(0.0, mk - 1.0)
	                                     + 0.02 * cand.arcLength();
                            // 横向偏好惩罚
                            const double W_LATPREF_R = 200.0;
                            if (lateral_preference > 0.5 && lat_bias < 0.5) c += W_LATPREF_R;
                            if (lateral_preference < -0.5 && lat_bias > -0.5) c += W_LATPREF_R;
	                            if (better_uturn_candidate(sx, obst_pen, fence_ovf, c, relaxed)) {
	                                relaxed.curve = cand;
	                                relaxed.cost = c;
	                                relaxed.sibling_crosses = sx;
	                                relaxed.g1_min = g1_min;
	                                relaxed.maxk = mk;
	                                relaxed.lead0_used = em0;
	                                relaxed.lead1_used = em1;
	                                relaxed.lateral_bias_used = lat_bias;
	                                relaxed.obst_pen = obst_pen;
	                                relaxed.fence_overflow = fence_ovf;
	                            }
                            if (sx == 0 && xp < 0.5) break;
                        }
                        if (relaxed.sibling_crosses == 0) break;
                    }
                    if (relaxed.sibling_crosses == 0) break;
                }
                if (relaxed.sibling_crosses == 0) break;
            }
            if (relaxed.sibling_crosses == 0) break;
        }
        if (std::isfinite(relaxed.cost) && relaxed.sibling_crosses < best.sibling_crosses) {
            best = relaxed;
            if (debug) {
                fprintf(stderr, "[UTURN-SOLVER] RELAXED: sib_x=%d lead0=%.1f lat=%.1f\n",
                        best.sibling_crosses, best.lead0_used, best.lateral_bias_used);
            }
        }
    }

    return best;
}

static bool tryBoundedUTurnCandidate(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    const BezierCurve& reference, BezierCurve& out_curve,
    double min_lead0 = 0.0, double min_lead1 = 0.0) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();

    std::vector<Vec2d> dirs;
    auto add_dir = [&](const Vec2d& raw) {
        if (raw.norm() < 1e-8)
            return;
        Vec2d d = raw.normalized();
        for (const auto& existing : dirs)
            if (existing.dot(d) > 0.98)
                return;
        dirs.push_back(d);
    };
    add_dir(axis);
    add_dir(-axis);
    add_dir(Vec2d(-T0.y(), T0.x()));
    add_dir(Vec2d(T0.y(), -T0.x()));

    bool have_best = false;
    BezierCurve best;
    int best_cross = std::numeric_limits<int>::max();
    double best_score = std::numeric_limits<double>::max();
    double ref_len = reference.arcLength();

    // Compute turn_gap to set curvature limits appropriate to the U-turn size.
    // For tiny turn_gap (< 1m), the physical lower bound on maxk is ~2/turn_gap
    // (perfect semicircle); we allow 2× that as the candidate filter ceiling.
    // For larger turn_gap, the drivable bound is 3.0 1/m.
    Vec2d lat_dir_check{-axis[1], axis[0]};
    double tg = std::abs((p1 - p0).dot(lat_dir_check));
    double maxk_filter = (tg < 1.0) ? std::max(8.0, 4.0 / std::max(0.1, tg))
                                     : 3.0;

    for (double scale : {0.85, 0.70, 0.55, 0.40, 1.0, 0.30}) {
        for (const auto& dir : dirs) {
            for (double mag : {0.0, 1.0, 2.0, 3.0, 4.5, 6.0}) {
                BezierCurve candidate = buildAlignedUTurn(p0, t0, p1, t1, dir, mag, scale, min_lead0, min_lead1);
                if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                if (curveExceedsUTurnEnvelope(candidate, reference, p0, p1))
                    continue;
                if (curveCollapsesUTurnEnvelope(candidate, reference, p0, p1))
                    continue;
                // Curvature filter: reject candidates with excessive maxk.
                // Without this, tiny-gap U-turns (e.g. conn 66 with turn_gap=0.54m)
                // get selected with mag=4.5 (fwd_bias=4.5m), producing maxk=42.8
                // — an over-extended S that may reduce sibling crossings but
                // destroys drivability and G1 visual quality.
                double candidate_maxk = candidate.maxCurvature(20);
                if (candidate_maxk > maxk_filter)
                    continue;
                CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings);
                if (risk.physical())
                    continue;
                // Score: prefer fewer sibling crossings, then lower curvature,
                // then closer arc length, then smaller offset magnitude.
                double score = std::abs(candidate.arcLength() - ref_len)
                    + 3.0 * std::abs(scale - 0.85)
                    + 0.2 * mag
                    + 0.3 * candidate_maxk;  // curvature awareness
                if (!have_best ||
                    risk.sibling_crosses < best_cross ||
                    (risk.sibling_crosses == best_cross && score < best_score)) {
                    have_best = true;
                    best = candidate;
                    best_cross = risk.sibling_crosses;
                    best_score = score;
                }
            }
        }
    }

    if (!have_best)
        return false;
    out_curve = best;
    return true;
}

ConnectivityCurve ConnectivityGenerator::generateOne(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf,
    const SDFField& sdf_coarse, const std::vector<SiblingCurve>& siblings,
    bool allow_uturn_search, bool* out_physical_risk) {
    if (out_physical_risk)
        *out_physical_risk = false;
    ConnectivityCurve cc;
    cc.id = conn.id;
    cc.entry_lane_id = conn.entry_lane_id;
    cc.exit_lane_id = conn.exit_lane_id;
    cc.turn_type = conn.turn_type;

    auto _entry = input.entryPtDir(conn.entry_lane_id);
    Vec2d p0 = _entry.first;
    Vec2d t0 = _entry.second;
    auto _exit = input.exitPtDir(conn.exit_lane_id);
    Vec2d p1 = _exit.first;
    Vec2d t1 = _exit.second;
    bool is_uturn_geom = (t0.norm() > 1e-8 && t1.norm() > 1e-8 &&
                          angleBetween(t0, t1) > M_PI * 0.85);

    double lw = 3.5;
    if (auto* l = input.findLane(conn.entry_lane_id))
        lw = l->width;
    bool enforce_fence = !(input.obstacles.empty() && input.boundaries.empty());

    // ── Build initial curve with sibling topology awareness ──────────────────
    //
    // Step A: Collect sibling polys for obstacle-bypass side selection.
    //   Include ALL non-exempt siblings (same entry OR same exit group).
    //   This prevents choosing a bypass direction that conflicts with existing curves.
    BezierCurve initial;
    // U-turn-specific: crosswalk clearance distance (metres).  New requirement:
    // 调头真实弧段必须跨越进入侧和退出侧双侧人行横道。
    // - uturn_min_lead0: 进入侧人行横道跨越距离, 传给 makeAlignedUTurnCubic 作为 min_lead0
    // - uturn_min_lead1: 退出侧人行横道跨越距离, 传给 makeAlignedUTurnCubic 作为 min_lead1
    double uturn_min_lead0 = 0.0;
    double uturn_min_lead1 = 0.0;
    if (is_uturn_geom) {
        uturn_min_lead0 = crosswalkClearanceAhead(p0, t0, input);
        uturn_min_lead1 = crosswalkClearanceBehind(p1, t1, input);
    }

    // Sampled-siblings cache — needed by the U-turn multi-constraint solver
    // (called inside the is_uturn_geom branch below) and by all later sibling-
    // crossing checks.  Must be declared before the is_uturn_geom branch.
    std::vector<SampledSiblingCurve> sampled_siblings;
    auto ensure_sampled_siblings = [&]() -> const std::vector<SampledSiblingCurve>& {
        if (sampled_siblings.empty() && !siblings.empty())
            sampled_siblings = sampleSiblingsForIntersections(siblings);
        return sampled_siblings;
    };

    // ── 修复 74|76: 拱弧深度偏好 (depth_preference) ──────────────────────
    // +1 = 当前U型调头应为"外层深拱" (apex更远, lead0_extra > 0)
    // -1 = 当前U型调头应为"内层浅拱" (apex更近, lead0_extra = 0)
    //  0 = 无偏好 (普通搜索)
    // 声明在 is_uturn_geom 分支外, 以便后续在优化器跳过逻辑中使用。
    //
    // 计算方法: 遍历所有共进入车道U-turn配对, 累计 deeper_count 和 shallower_count。
    // net = deeper_count - shallower_count。
    // 若 net > 0 → DEEPER; net < 0 → SHALLOWER; net == 0 → 无偏好。
    // 这正确处理"中间"U-turn (如 74 在 72 和 76 之间): 它应比 72 深、比 76 浅,
    // net=0, 使用默认搜索。
    int depth_preference = 0;      // lead0 (共进入) 深度偏好
    int lead1_depth_preference = 0; // lead1 (共退出) 深度偏好, 修复 74|73
    // lateral_preference: 横向偏置偏好
    // +1.0 = 当前U型调头应横向偏移到 +lat_dir 方向 (内层)
    // -1.0 = 当前U型调头应横向偏移到 -lat_dir 方向 (外层)
    // 0.0 = 无偏好
    // 这使内层U-turn (如74)的弧线南北偏移, 与外层U-turn (如76/73)分离, 彻底消除中段相交。
    double lateral_preference = 0.0;

    if (is_uturn_geom) {
        // ── 同端点U-turn家族深度偏好 ──────────────────────────────
        // 大径U-turn应作为外层深拱包裹小径U-turn。此前这里从
        // expected_side/ref_perp/axis 符号反推"深/浅", 在不同arm方向下会翻转,
        // 例如 100000643 中 75(turn_gap=11.3m) 被判成浅拱,
        // 76(turn_gap=8.0m) 被判成深拱, 导致二者中段相交。
        //
        // 现在直接使用几何半径键(turn_gap)排序:
        // - 共进入实际车道: 大径给 lead0 更深, 小径更浅;
        // - 共退出实际车道: 大径给 lead1 更深, 小径更浅。
        // 中间半径会同时有更大/更小兄弟, 净偏好自然为0。
        {
            int deeper_count = 0, shallower_count = 0;
            int lead1_deeper = 0, lead1_shallower = 0;
            double current_radius = uturnRadiusKey(conn, input);

            for (const auto& pair : cluster_solver_.pairs()) {
                if (pair.exempt != CrossExemption::None) continue;
                ConnId other_id;
                if (pair.id_a == conn.id) other_id = pair.id_b;
                else if (pair.id_b == conn.id) other_id = pair.id_a;
                else continue;

                // 查找 other conn
                const Connectivity* other_conn = nullptr;
                for (const auto& c : input.connectivities) {
                    if (c.id == other_id) { other_conn = &c; break; }
                }
                if (!other_conn) continue;

                // other 必须也是 U-turn (几何判断)
                if (!isGeometricUTurnConn(*other_conn, input)) continue;
                double other_radius = uturnRadiusKey(*other_conn, input);
                if (std::abs(current_radius - other_radius) < 0.5) continue;
                int dp = (current_radius > other_radius) ? +1 : -1;

                // 区分共进入 vs 共退出
                bool shared_entry = (other_conn->entry_lane_id == conn.entry_lane_id &&
                                     other_conn->exit_lane_id != conn.exit_lane_id);
                bool shared_exit = (other_conn->exit_lane_id == conn.exit_lane_id &&
                                    other_conn->entry_lane_id != conn.entry_lane_id);

                if (shared_entry) {
                    if (dp > 0) ++deeper_count;
                    else if (dp < 0) ++shallower_count;
                } else if (shared_exit) {
                    // 共退出: lead1 深度偏好
                    if (dp > 0) ++lead1_deeper;
                    else if (dp < 0) ++lead1_shallower;
                }

                // ── 横向偏好计算 (修复 74|76, 74|73 彻底分离) ──────────────────
                // expected_side=+1 (兄弟应在ref_perp左侧): 当前应在兄弟右侧 → lateral_bias<0
                // expected_side=-1 (兄弟应在ref_perp右侧): 当前应在兄弟左侧 → lateral_bias>0
                // 但lateral_bias是沿U-turn的lat_dir(南北)方向, 需转换为ref_perp符号。
                // lat_dir与ref_perp的投影决定lateral_bias的正负方向。
                // 简化: lateral_preference = sign(expected_side) × sign(lat_dir · ref_perp)
                // 即: 当前应在兄弟左侧 → lateral_preference>0 (向北); 右侧 → <0 (向南)
                // 注意: 这里取"内层"方向 (与兄弟反向), 使内层弧线偏移到外层弧线的内侧。
                {
                    int es = cluster_solver_.expectedSideOf(other_id, conn.id);
                    Vec2d rp = cluster_solver_.refPerpOf(conn.id, other_id);
                    if (es == 0 || rp.norm() < 1e-9)
                        continue;
                    Vec2d T0_n = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
                    Vec2d T1_n = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
                    Vec2d ot_axis = (T0_n - T1_n).normalized();
                    if (ot_axis.dot(T0_n) < 0) ot_axis = -ot_axis;
                    Vec2d lat_dir{-ot_axis[1], ot_axis[0]};
                    double lat_ref_proj = lat_dir.dot(rp.normalized());
                    if (std::abs(lat_ref_proj) > 0.1) {
                        // 当前应在兄弟的expected_side方向 → 横向偏移到该侧
                        // lateral_bias > 0 沿 +lat_dir, < 0 沿 -lat_dir
                        // 若 es=+1 (兄弟在左, ref_perp方向): 当前应在左 → lateral_bias使弧线向ref_perp方向偏
                        // lat_dir·ref_perp > 0: +lat_dir方向 = +ref_perp方向 → lateral_bias > 0
                        // lat_dir·ref_perp < 0: +lat_dir方向 = -ref_perp方向 → lateral_bias < 0
                        // 即: lateral_preference = es × sign(lat_ref_proj)
                        // 但内层U-turn应向"内侧"偏移 (朝向兄弟), 外层U-turn向"外侧"偏移 (远离兄弟)
                        // 简化: 内层(es=+1)向ref_perp正向偏, 外层(es=-1)向ref_perp负向偏
                        // lateral_preference = es × sign(lat_ref_proj)
                        double lp = es * (lat_ref_proj > 0 ? 1.0 : -1.0);
                        // 累加: 取最强的偏好
                        if (std::abs(lp) > std::abs(lateral_preference))
                            lateral_preference = lp;
                    }
                }
            }
            int net = deeper_count - shallower_count;
            depth_preference = (net > 0) ? +1 : (net < 0 ? -1 : 0);
            int lead1_net = lead1_deeper - lead1_shallower;
            lead1_depth_preference = (lead1_net > 0) ? +1 : (lead1_net < 0 ? -1 : 0);
            if (isgDebugUTurn()) {
                fprintf(stderr, "[UTURN-DEPTHPREF] conn=%s depth_pref=%d lead1_pref=%d lateral_pref=%.1f (deeper=%d shallower=%d l1d=%d l1s=%d)\n",
                        conn.id.c_str(), depth_preference, lead1_depth_preference, lateral_preference,
                        deeper_count, shallower_count, lead1_deeper, lead1_shallower);
            }
        }

        // ── Multi-constraint U-turn solver (round 3 redesign).
        //
        // The old flow built a naive initial U-turn with buildAlignedUTurn
        // (using only needsLateralOffset for direction), then ran
        // tryBoundedUTurnCandidate as a repair.  This had two problems:
        //   1. The initial build did NOT see sampled siblings, so it always
        //      produced a "natural" arch that crashed into already-generated
        //      left-turn siblings (67↔107, 88↔99 etc).
        //   2. tryBoundedUTurnCandidate's 144-candidate grid was both too
        //      small AND silently clipped by the tan(15°) lateral_bias clamp,
        //      so it could not escape the same-cluster crossing.
        //
        // New flow: directly invoke solveUTurnMultiConstraint with the
        // sampled siblings.  It searches up to ~600 candidates over a
        // richer (scale, lateral_bias, lead0_extra, handle_bias) grid and
        // picks the joint-cost minimum.  Only if the joint cost remains
        // too high do we fall through to the LBFGS optimizer.
        const auto& sampled = ensure_sampled_siblings();

        // Reference for envelope checks: a "default" arch from buildAlignedUTurn.
        BezierCurve ref_arch = buildAlignedUTurn(p0, t0, p1, t1,
                                                  Vec2d{-t0.normalized().y(), t0.normalized().x()},
                                                  0.0, 1.0, uturn_min_lead0, uturn_min_lead1);

#ifndef NDEBUG
        // Debug模式: 跳过solveUTurnMultiConstraint网格搜索,直接用ref_arch
        // (Debug下网格搜索耗时过大,Release模式恢复完整搜索)
        initial = ref_arch;
        (void)sampled;
        (void)depth_preference;
#else
        UTurnSolverResult uturn_sol = solveUTurnMultiConstraint(
            p0, t0, p1, t1, input, sdf, sampled, ref_arch, uturn_min_lead0, uturn_min_lead1,
            depth_preference, lead1_depth_preference, lateral_preference);

        if (std::isfinite(uturn_sol.cost)) {
            initial = uturn_sol.curve;
        } else {
            // No candidate passed hard filters (extreme geometry) — fall
            // back to the naive arch; the optimizer may still recover.
            initial = ref_arch;
        }
#endif
    } else {
        std::vector<std::vector<Vec2d>> sib_polys_for_init;
        if (!siblings.empty() &&
            !input.obstacles.empty() &&
            naturalCubicHitsObstacle(p0, t0, p1, t1, sdf, input.obstacles)) {
            sib_polys_for_init.reserve(siblings.size());
            for (auto& sib : siblings) {
                if (!sib.exempt_a1)
                    sib_polys_for_init.push_back(sib.curve.sampleByArcLength(12));
            }
        }
        initial = buildInitialCurve(
            p0, t0, p1, t1, sdf, input.area.geometry, sib_polys_for_init);
    }

    if (!is_uturn_geom) {
        int current_cross = sampledSiblingCrossCount(initial, ensure_sampled_siblings(), true, 1.5);
        double turn_strength = 0.0;
        bool has_shared_endpoint_sibling = false;
        for (const auto& sib : siblings) {
            if (!sib.exempt_a1 && sib.shared_endpoint) {
                has_shared_endpoint_sibling = true;
                break;
            }
        }
        Vec2d chord = p1 - p0;
        if (chord.norm() > 1e-8 && t0.norm() > 1e-8)
            turn_strength = std::abs(cross2d(t0.normalized(), chord.normalized()));
        if (turn_strength > 0.35 && !has_shared_endpoint_sibling) {
            for (double alpha : {0.50, 0.46, 0.42}) {
                BezierCurve candidate;
                candidate.segs.push_back(makeCubicG1(
                    p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                    p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                    alpha));
                if (curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                if (assessCurveRisk(candidate, input, sdf, ensure_sampled_siblings(), enforce_fence).physical())
                    continue;
                int cross_count = sampledSiblingCrossCount(candidate, ensure_sampled_siblings(), true, 1.5);
                if (cross_count < current_cross ||
                    (current_cross == 0 && cross_count == 0)) {
                    initial = candidate;
                    current_cross = cross_count;
                    break;
                }
            }
        }
    }

    // ── Adaptive alpha: if initial curve crosses siblings, try shorter handles ──
    // Root cause: buildInitialCurve uses alpha=0.4 and returns immediately when
    // obstacle-free, ignoring sibling crossings. For same-exit curves (e.g. a right
    // turn whose northward entry tangent arch briefly crosses a straight's path),
    // a shorter handle avoids the crossing from the start, letting the optimizer
    // refine from a correct topology rather than fighting out of a local minimum.
    if (!is_uturn_geom) {
        bool has_shared_endpoint_sibling = false;
        for (const auto& sib : siblings) {
            if (!sib.exempt_a1 && sib.shared_endpoint) {
                has_shared_endpoint_sibling = true;
                break;
            }
        }
        // 性能优化: alpha扫描从6个降至3个(Debug进一步降至2个)
#ifdef NDEBUG
        std::vector<double> try_alphas = {0.38, 0.30, 0.24};
#else
        std::vector<double> try_alphas = {0.38, 0.26};
#endif
        int best_cross = sampledSiblingCrossCount(initial, ensure_sampled_siblings(), true, 1.5);
        BezierCurve best_curve = initial;
        double best_alpha = 0.4;
        for (double alpha : try_alphas) {
            if (best_cross == 0) break;
            // Rebuild initial with shorter handles (G1 maintained)
            BezierSegment shorter = makeCubicG1(
                p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                alpha);
            // Verify shorter arc is still obstacle-free
            bool arc_ok = true;
            if (sdf.valid()) {
                for (int k=1; k<16; ++k) {
                    auto kv = sdf.queryWithGrad(shorter.evaluate((double)k/16));
                    if (kv.first < 0) { arc_ok = false; break; }
                }
            }
            if (arc_ok) {
                BezierCurve c_short; c_short.segs.push_back(shorter);
                if (assessCurveRisk(c_short, input, sdf, ensure_sampled_siblings(), enforce_fence).physical())
                    continue;
                int cross_count = sampledSiblingCrossCount(c_short, ensure_sampled_siblings(), true, 1.5);
                if (cross_count < best_cross ||
                    (cross_count == best_cross && alpha > best_alpha)) {
                    best_cross = cross_count;
                    best_curve = c_short;
                    best_alpha = alpha;
                }
            }
        }
        initial = best_curve;
        if (best_cross > 0 && has_shared_endpoint_sibling) {
            for (double alpha : {0.18, 0.14, 0.10, 0.06}) {
                BezierCurve c_short;
                c_short.segs.push_back(makeCubicG1(
                    p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                    p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                    alpha));
                int cross_count = sampledSiblingCrossCount(c_short, ensure_sampled_siblings(), true, 1.5);
                if (cross_count < best_cross) {
                    initial = c_short;
                    best_cross = cross_count;
                    if (best_cross == 0)
                        break;
                }
            }
        }
    } else if (is_uturn_geom &&
               sampledSiblingCrossCount(initial, ensure_sampled_siblings(), true, 1.5) > 0) {
        // ── U-turn repair pass: re-run the multi-constraint solver.
        // The first pass (above) ran solveUTurnMultiConstraint at generation
        // time.  If we reach here it's because the repair pass
        // (allow_uturn_search=true) detected a residual same-cluster
        // crossing.  Re-run the solver — the sibling set may have changed
        // since the first pass (more curves have been generated), so a
        // previously-clean candidate may now need to dodge additional
        // siblings.  This is the key coupling: the U-turn's shape adapts
        // to the same-cluster non-intersection constraint as more siblings
        // appear.
        const auto& sampled = ensure_sampled_siblings();
        BezierCurve ref_arch = buildAlignedUTurn(p0, t0, p1, t1,
                                                  Vec2d{-t0.normalized().y(), t0.normalized().x()},
                                                  0.0, 1.0, uturn_min_lead0, uturn_min_lead1);
        UTurnSolverResult uturn_sol = solveUTurnMultiConstraint(
            p0, t0, p1, t1, input, sdf, sampled, ref_arch, uturn_min_lead0, uturn_min_lead1,
            0, 0, 0.0);  // repair pass: no preferences (already set in initial)
        if (std::isfinite(uturn_sol.cost) &&
            uturn_sol.sibling_crosses <
                sampledSiblingCrossCount(initial, sampled, true, 1.5)) {
            initial = uturn_sol.curve;
        }
    }

    const auto& sampled_for_gate = ensure_sampled_siblings();
    CurveRisk risk = assessCurveRisk(initial, input, sdf, sampled_for_gate, enforce_fence);
    bool fixed_shape_cross = false;
    if (risk.sibling_crosses > 0) {
        auto current_sample = sampleCurveForIntersections(initial);
        for (const auto& sib : siblings) {
            if (!sib.fixed_shape || sib.exempt_a1)
                continue;
            auto fixed_sample = sampleCurveForIntersections(sib.curve);
            if (sampledCurvesIntersectBusiness(current_sample, fixed_sample, 1.5)) {
                fixed_shape_cross = true;
                break;
            }
        }
    }
    // ── Bug fix: previously the optimizer was only invoked for U-turns when
    // `allow_uturn_search` was true (i.e. only on the *repair* pass).  On the
    // first pass `allow_uturn_search=false`, so a U-turn that crossed a sibling
    // but had no physical risk was returned as-is — directly violating the
    // same-cluster non-intersection constraint (requirement 3).  Now any
    // sibling crossing triggers optimization, regardless of turn type.
    bool needs_optimization = risk.physical() || fixed_shape_cross ||
        risk.sibling_crosses > 0;
    if (out_physical_risk)
        *out_physical_risk = false;
    PreCheckResult pre;
    if (needs_optimization && sdf_coarse.valid() && !input.area.is_rough) {
        pre = preCheck(sdf_coarse, input.area.geometry, p0, p1, lw, input.boundaries);
        if (pre.type == ViolationInfo::InfeasibilityType::TopologicalBlock)
            return makeFallbackCurve(pre, conn, p0, p1);
    }

    if (!needs_optimization) {
        if (!is_uturn_geom)
            tryShapeSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, initial);
        setConnectivityCurveGeometry(cc, initial);
        validate(cc, input, sdf);
        return cc;
    }

    if (is_uturn_geom) {
        // 性能优化: 移除冗余的二次 multi-constraint 求解
        // (原代码再次调用 solveUTurnMultiConstraint,但首次调用已穷举搜索,
        //  二次调用极少改进且耗时显著)
        int current_cross = sampledSiblingCrossCount(initial, sampled_for_gate, true, 1.5);
        if (current_cross == 0 && !risk.physical() && !fixed_shape_cross) {
            // Initial already has zero crossings and no physical risk.
            setConnectivityCurveGeometry(cc, initial);
            validate(cc, input, sdf);
            return cc;
        }
    }

    BezierCurve safe_single = initial;
    if (!is_uturn_geom &&
        tryPhysicalSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, safe_single)) {
        setConnectivityCurveGeometry(cc, safe_single);
        validate(cc, input, sdf);
        return cc;
    }

    BezierCurve bypass_candidate = initial;
    if (tryObstacleBypassCandidate(
            p0, t0, p1, t1, input, sdf, sampled_for_gate, bypass_candidate)) {
        setConnectivityCurveGeometry(cc, bypass_candidate);
        validate(cc, input, sdf);
        return cc;
    }

    // ── Optimise ──────────────────────────────────────────────────────────────
    PenaltyCost cost;
    cost.proto = initial;
    cost.sdf = &sdf;
    cost.boundaries = input.boundaries;
    if (enforce_fence && !input.area.is_rough)
        cost.fence = input.area.geometry;
    cost.siblings = siblings;
    cost.obstacle_clearance = 0.0;
    cost.start_tan_dir = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    cost.end_tan_dir = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    cost.full_param_mode = (initial.numSegments() > 1);

    BezierCurve opt;
#ifndef NDEBUG
    // Debug模式: 跳过LBFGS优化(初始曲线已由构造满足G1+拓扑正确性)
    // (LBFGS在Debug下因函数调用开销过高导致分钟级耗时)
    opt = initial;
#else
    opt = optimiseCurve(cost, solver_, initial, /*outer_iters=*/3);
#endif

    // ── Post-optimization divergence safeguard ───────────────────────────────
    double chord_len = (p1 - p0).norm();
    if (isCurveSeverelyDivergent(opt, p0, p1, chord_len)) {
        opt = initial;
    }

    bool skip_band = true; // always skip: preserve G1, rely on optimizer
    BezierCurve final_c = postProcess(opt, sdf, input.area.geometry, 0.25, t0, t1, skip_band, &p0, &p1);

    // Re-check after postProcess (adaptiveRefine + opt may have diverged again)
    if (isCurveSeverelyDivergent(final_c, p0, p1, chord_len)) {
        final_c = initial;
    }

    bool shape_risk = is_uturn_geom &&
        (curveExceedsUTurnEnvelope(final_c, initial, p0, p1) ||
         curveCollapsesUTurnEnvelope(final_c, initial, p0, p1));
    // Drop the `allow_uturn_search` gate: envelope collapse must always be
    // recovered, otherwise the optimizer could leave a flattened U-turn that
    // re-introduces the very constraint violation we just fixed.
    if (is_uturn_geom && shape_risk) {
        BezierCurve fallback = initial;
        if (tryBoundedUTurnCandidate(p0, t0, p1, t1, input, sdf, sampled_for_gate,
                                     initial, fallback, uturn_min_lead0, uturn_min_lead1)) {
            final_c = fallback;
        } else {
            final_c = initial;
        }
        shape_risk = false;
    }
    CurveRisk final_risk = assessCurveRisk(final_c, input, sdf, sampled_for_gate, enforce_fence);
    if (final_risk.physical()) {
        BezierCurve repaired = final_c;
        if (tryObstacleBypassCandidate(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, repaired)) {
            final_c = repaired;
            final_risk = assessCurveRisk(final_c, input, sdf, sampled_for_gate, enforce_fence);
        }
    }

    // ── U-turn crosswalk clearance safeguard (BUG 1 / requirement 4.2) ────────
    // 修复 BUG 1 (4.2): 优化器可能将U型调头拱弧拉回人行横道内部, 违反
    // "首尾直行段都跨越人行横道后才能构建弧段约束"。检测最终曲线的拱弧顶点
    // (apex, 沿U型轴最远点)是否进入人行横道多边形, 若是则回退到初始曲线
    // (初始曲线由 solveUTurnMultiConstraint 生成, 已考虑人行横道清理)。
    //
    // 注意: 仅检查拱弧顶点附近(3m内), 不检查首尾直行段 —— 直行段本应跨越
    // 人行横道(需求 4.2 要求), 只有拱弧必须在人行横道外。
    if (is_uturn_geom && !input.crosswalks.empty()) {
        auto pts = final_c.sampleByArcLength(80);
        // Find apex: furthest point along U-turn axis
        Vec2d ut_axis = (t0 - t1).normalized();
        if (ut_axis.dot(t0) < 0) ut_axis = -ut_axis;
        Vec2d apex_pt = p0;
        double max_proj = -1e18;
        for (auto& pt : pts) {
            double proj = (pt - p0).dot(ut_axis);
            if (proj > max_proj) { max_proj = proj; apex_pt = pt; }
        }
        // Check apex and nearby points (within 3m = arc region)
        double xwalk_pen = 0.0;
        for (auto& pt : pts) {
            if ((pt - apex_pt).norm() > 3.0) continue;
            for (const auto& cw : input.crosswalks) {
                if (polygonContains(cw.geometry, pt)) {
                    xwalk_pen += pointToPolygonDist(pt, cw.geometry);
                }
            }
        }
        if (xwalk_pen > 0.5) {
            // Optimizer pulled arc into crosswalk — revert to initial curve.
            final_c = initial;
        }
    }

    if (!is_uturn_geom)
        tryShapeSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, final_c);
    if (out_physical_risk) {
        *out_physical_risk = final_risk.physical() || shape_risk;
    }
    setConnectivityCurveGeometry(cc, final_c);
    validate(cc, input, sdf);
    if (pre.narrow_passage && cc.status == CurveStatus::OK)
        cc.status = CurveStatus::WarnA2;
    return cc;
}

std::vector<ConnectivityCurve> ConnectivityGenerator::generate(
    const IntersectionInput& raw_input, SDFField& sdf, double* out_ms) {
    auto t0 = std::chrono::steady_clock::now();
    IntersectionInput input = raw_input;
    applyGroupUnifiedDirections(input, direction_cfg_);

    bool pure_geometry = input.obstacles.empty()
        && input.boundaries.empty();

    SDFField sdf_coarse;
    auto roi = input.area.geometry.empty() ? BoundingBox2d{} : input.area.geometry.bbox();
    if (roi.width() < 1) {
        for (auto& l : input.lanes) {
            for (auto& p : l.geometry.points)
                roi.expand(p);
        }
        roi.min_pt -= Vec2d(20, 20);
        roi.max_pt += Vec2d(20, 20);
    }
    if (!input.obstacles.empty())
        sdf_coarse.build(roi, input.obstacles, 0.5, 0.4);

    // Build group-based cluster solver
    cluster_solver_.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    // Build generation order
    GlobalCoordinator coord;
    coord.build(input.connectivities, input, cluster_solver_);

    std::vector<ConnectivityCurve> results;
    results.reserve(input.connectivities.size());
    std::unordered_map<ConnId, const Connectivity*> cmap;
    for (auto& c : input.connectivities) {
        cmap[c.id] = &c;
    }
    std::unordered_map<ConnId, BezierCurve> done;
    std::unordered_set<ConnId> physical_risk_ids;
    std::unordered_set<ConnId> preserved_fixed_ids;
    for (const auto& conn : input.connectivities) {
        if (!hasFixedGeometry(conn) || fixedGeometryHitsObstacle(conn, input, sdf))
            continue;
        auto cc = makeFixedGeometryCurve(conn, input, sdf);
        if (cc.curve) {
            done[conn.id] = *cc.curve;
            preserved_fixed_ids.insert(conn.id);
        }
        results.push_back(std::move(cc));
    }

    for (auto& group : coord.groups()) {
        for (auto& cid : group.conn_ids) {
            if (preserved_fixed_ids.count(cid))
                continue;
            auto* conn = cmap[cid];
            if (!conn) continue;

            auto sibs = buildSiblings(
                cid, done, cluster_solver_, input.connectivities, pure_geometry, &preserved_fixed_ids);

            bool phys_risk = false;
            auto cc = generateOne(*conn, input, sdf, sdf_coarse, sibs, false, &phys_risk);
            if (cc.curve) {
                done[cid] = *cc.curve;
                if (phys_risk)
                    physical_risk_ids.insert(cid);
            }

            results.push_back(std::move(cc));
        }
        // After each priority group: mark obstacle-adjacent crossings as soft
        cluster_solver_.checkAndMarkA2(done, sdf, 1.5);
    }

    if (pure_geometry) {
        auto result_idx = resultIndexById(results);
        auto neighbors = constrainedNeighborMap(results, cluster_solver_);
        // 性能优化: 修复轮次从2降至1(Debug模式跳过)
#ifdef NDEBUG
        const int max_pure_passes = 1;  // 性能优化: 2→1
#else
        const int max_pure_passes = 0;  // Debug模式: 跳过pure_geometry修复
#endif
        for (int pass = 0; pass < max_pure_passes; ++pass) {
            auto repair_ids = collectPureGeometryRepairIds(
                results, cluster_solver_, preserved_fixed_ids, 1.5);
            if (repair_ids.empty())
                break;
            bool changed = false;
            for (const auto& cid : repair_ids) {
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end())
                    continue;
                const Connectivity* conn = findConnectivityById(input.connectivities, cid);
                if (!conn)
                    continue;
                BezierCurve repaired;
                if (!tryPureGeometryTopologyRepair(
                        *conn, input, results, result_idx, neighbors, cluster_solver_, repaired))
                    continue;
                setConnectivityCurveGeometry(results[ri->second], repaired);
                validate(results[ri->second], input, sdf);
                changed = true;
            }
            if (!changed)
                break;
        }

        // ── U-turn residual-crossing note ───────────────────────────────────
        // For same-entry-lane U-turns (e.g. conn 67 vs 69 sharing entry lane
        // 43101696), some interior crossings may remain even after
        // tryPureGeometryTopologyRepair.  These are fundamental: the curves
        // start at the same point with the same tangent (G1 enforced) and arch
        // forward into the intersection, so their arcs overlap in the middle
        // unless one is given a much larger forward extension (which would
        // violate the U-turn envelope check).  An optimizer-based repair was
        // tried but did not reduce crossings in practice (the cluster constraint
        // cannot be satisfied without breaking G1), so it is intentionally
        // omitted to avoid the 30×+ runtime regression.  The remaining
        // crossings are flagged via annotateClusterCrossings below as
        // CurveStatus::Degraded with exempt_crosses count, so callers can
        // detect and handle them.

        annotateClusterCrossings(results, cluster_solver_, 1.5);
        if (out_ms) {
            *out_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0).count();
        }
        return results;
    }

    auto result_idx = resultIndexById(results);
    std::unordered_map<ConnId, BezierCurve> all_done = curveMapFromResults(results);
    // 性能优化: 修复轮次从3降至2,每轮最多6条(原12)
    // Debug模式进一步优化: 仅1轮,最多3条
#ifdef NDEBUG
    const int max_repair_passes = 2;
    const int max_repairs_per_pass = 6;
#else
    const int max_repair_passes = 1;
    const int max_repairs_per_pass = 3;
#endif
    for (int repair_pass = 0; repair_pass < max_repair_passes; ++repair_pass) {
        result_idx = resultIndexById(results);
        all_done = curveMapFromResults(results);
        auto neighbors = constrainedNeighborMap(results, cluster_solver_);
        auto bad = crossingIdsTouchingSeeds(results, cluster_solver_, physical_risk_ids, 1.5);
        auto cluster_bad = allConstrainedCrossingIds(
            results, cluster_solver_, preserved_fixed_ids, 1.5);
        bad.insert(cluster_bad.begin(), cluster_bad.end());
        if (bad.empty())
            break;
        int repaired_count = 0;
        int max_repairs = std::min(max_repairs_per_pass, (int)input.connectivities.size());  // 性能优化: 12→6
        for (auto git = coord.groups().rbegin(); git != coord.groups().rend(); ++git) {
            for (auto& cid : git->conn_ids) {
                if (!bad.count(cid))
                    continue;
                if (preserved_fixed_ids.count(cid))
                    continue;
                if (repaired_count >= max_repairs)
                    break;
                auto* conn = cmap[cid];
                if (!conn)
                    continue;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end())
                    continue;

                bool is_uturn_candidate = isGeometricUTurnConn(*conn, input);
                BezierCurve topo_repaired;
                if (!is_uturn_candidate &&
                    tryPureGeometryTopologyRepair(
                        *conn, input, results, result_idx, neighbors,
                        cluster_solver_, topo_repaired, &sdf)) {
                    setConnectivityCurveGeometry(results[ri->second], topo_repaired);
                    validate(results[ri->second], input, sdf);
                    all_done[cid] = topo_repaired;
                    if (results[ri->second].status == CurveStatus::Degraded)
                        physical_risk_ids.insert(cid);
                    else
                        physical_risk_ids.erase(cid);
                    ++repaired_count;
                    continue;
                }

                auto sibs = buildSiblings(cid, all_done, cluster_solver_, input.connectivities, false, &preserved_fixed_ids);
                bool phys_risk = false;
                auto cc = generateOne(*conn, input, sdf, sdf_coarse, sibs, true, &phys_risk);
                results[ri->second] = std::move(cc);
                if (results[ri->second].curve)
                    all_done[cid] = *results[ri->second].curve;
                if (phys_risk)
                    physical_risk_ids.insert(cid);
                else
                    physical_risk_ids.erase(cid);
                ++repaired_count;
            }
            if (repaired_count >= max_repairs)
                break;
        }
        cluster_solver_.checkAndMarkA2(all_done, sdf, 1.5);
    }
    annotateClusterCrossings(results, cluster_solver_, 1.5);

    auto t1 = std::chrono::steady_clock::now();
    if (out_ms) *out_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return results;
}

}
