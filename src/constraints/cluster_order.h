#pragma once

#include "types.h"
#include <unordered_map>

namespace isg {

enum class CrossExemption {
    None, StructuralCross, ObstacleCross
};

struct CurvePair {
    ConnId id_a, id_b;
    CrossExemption exempt = CrossExemption::None; //交叉豁免
    double exempt_zone_radius = 0.0; //豁免区半径
    int expected_side = 0;   // 预期方位, +1: id_a在左侧; -1: id_a在右侧; 0: unknown
    Vec2d ref_perp{0, 0}; // pair-specific: left-normal of mean_entry→mean_exit
    // shared_endpoint: both curves go to the SAME exit lane (same entry group).
    // Ordering is maintained but the cluster penalty uses a wider skip zone near
    // the convergence point so the constraint is not applied at the shared endpoint.
    bool shared_endpoint = false;
};

class SDFField;

// ─────────────────────────────────────────────────────────────────────────────
// ClusterOrderSolver — Topology-first design
//
// Sort strategy:
//   Entry cluster: DESCENDING angle(exit_pt - mean_entry_pt)
//     tie-break by exit lane lateral position within exit group (geometric)
//   Exit cluster: DESCENDING projection of entry_pt onto exit arm left-normal
//     tie-break by exit lane lateral position within exit group (geometric)
//
// Structural cross detection: topological rank × lateral inversion
//   For same-entry-group pair (A,B):
//     entry_cluster_rank(A) < entry_cluster_rank(B) [A is LEFT in entry cluster]
//     but exit_lat(A) < exit_lat(B) [A exits MORE RIGHT]
//     → (rank_A-rank_B)*(exit_lat_A-exit_lat_B) > 0 → StructuralCross
//   For same-exit-group cross-arm pair (A,B):
//     exit_cluster_rank(A) < exit_cluster_rank(B) [A is LEFT in exit cluster]
//     but entry_lat_A < entry_lat_B [A approaches from RIGHT]
//     → (rank_A-rank_B)*(entry_lat_A-entry_lat_B) > 0 → StructuralCross
// ─────────────────────────────────────────────────────────────────────────────
class ClusterOrderSolver {
public:
    void build(const std::vector<Connectivity>&, const std::vector<Lane>&, const std::vector<LaneGroup>&);

    void markObstacleExempt(CurvePair&, const Vec2d&, const SDFField&, double r = 1.5);

    const std::vector<CurvePair>& pairs() const { return pairs_; }

    // Returns true when the pair (a,b) shares the same exit lane (same entry group).
    // Used to apply wider skip zone in evalCluster.
    bool isSharedEndpoint(const ConnId& a, const ConnId& b) const {
        for (auto& p : pairs_) {
            if ((p.id_a==a && p.id_b==b) || (p.id_a==b && p.id_b==a))
                return p.shared_endpoint;
        }
        return false;
    }

    CrossExemption exemptionOf(const ConnId&, const ConnId&) const;

    bool pairExists(const ConnId& a, const ConnId& b) const;

    int expectedSideOf(const ConnId &, const ConnId &) const;

    Vec2d refPerpOf(const ConnId &, const ConnId &) const;

    void checkAndMarkA2(const std::unordered_map<ConnId, BezierCurve>&, const SDFField&, double r = 1.5);

    const std::unordered_map<LaneGroupId, std::vector<ConnId>>& entryGroupOrder() const {
        return entry_group_order_;
    }

    const std::unordered_map<LaneGroupId, std::vector<ConnId>>& exitGroupOrder() const {
        return exit_group_order_;
    }

private:
    void addPairsFromSortedCluster(
            const std::vector<ConnId>& cids, const std::vector<Connectivity>& conns, const std::vector<Lane> &);

    static bool hasPair(const std::vector<CurvePair>& pairs, const ConnId& a, const ConnId& b);

    void detectTopologicalInversions(const std::vector<Connectivity> &);

private:
    std::vector<CurvePair> pairs_;
    std::unordered_map<LaneGroupId, std::vector<ConnId>> entry_group_order_;
    std::unordered_map<LaneGroupId, std::vector<ConnId>> exit_group_order_;

    // U-turn detection via t0·t1 dot product (not group ID comparison
    // which fails when enterGroupId ≠ exitGroupId even for physical U-turns).
    // Populated in build() from actual lane tangents.
    std::unordered_map<ConnId, bool> is_uturn_;

    // exit_lat_in_entry_ref sign: +1 = exits LEFT (west for south arm), -1 = exits RIGHT.
    // Used to detect opposite-side structural crosses (left turn vs right turn always cross).
    std::unordered_map<ConnId, double> exit_lat_sign_;

    // entry lateral position within entry group frame.
    // Used for entry-exit inversion detection:
    //   (entry_lat_A - entry_lat_B) × (exit_lat_A - exit_lat_B) < 0 → StructuralCross
    // Catches same-side exits with inverted lane ordering
    // (e.g. outer right turn from inner entry lane = must cross inner right turn from outer lane).
    std::unordered_map<ConnId, double> entry_lat_in_entry_ref_;

    // Lateral positions used by topological inversion detector
    // entry_cluster: exit_lat_in_entry_ref[cid] = exit_pt projected onto entry arm perp
    // exit_cluster:  entry_lat_in_exit_ref[cid]  = entry_pt projected onto exit arm left-normal
    std::unordered_map<ConnId, double> exit_lat_in_entry_ref_;
    std::unordered_map<ConnId, double> entry_lat_in_exit_ref_;

    // exit lane's lateral position in its own exit group frame.
    // = (exitEndpoint - mean_exit_pt) · exit_arm_left_normal
    // Distinct from entry_lat_in_exit_ref_ (which is always consistent with the
    // exit_cluster_rank_, making inversion detection a no-op for same_exit pairs).
    std::unordered_map<ConnId, double> exit_lat_in_exit_group_ref_;
    // Cluster ranks (lower index = more LEFT)
    std::unordered_map<ConnId, int> entry_cluster_rank_;
    std::unordered_map<ConnId, int> exit_cluster_rank_;
};

}