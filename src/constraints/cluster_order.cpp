#include "cluster_order.h"
#include "optimizer/sdf_field.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace isg {

// ─── Lane geometry helpers ────────────────────────────────────────────────────
static const Connectivity *findConn(const std::vector<Connectivity> &cs, const ConnId &id) {
    for (auto &c:cs) if (c.id == id)return &c;
    return nullptr;
}

static const Lane *findLane(const std::vector<Lane> &ls, const LaneId &id) {
    for (auto &l:ls) if (l.id == id)return &l;
    return nullptr;
}

static Vec2d entryEndpoint(const Lane *l) {
    if (!l || l->geometry.points.empty())return Vec2d(0, 0);
    return l->geometry.points.back();
}

static Vec2d exitEndpoint(const Lane *l) {
    if (!l || l->geometry.points.empty())return Vec2d(0, 0);
    return l->geometry.points.front();
}

static Vec2d entryTangent(const Lane *l) {
    if (!l || l->geometry.points.size() < 2)return Vec2d(0, 1);
    auto &p = l->geometry.points;
    Vec2d d = p.back() - p[p.size() - 2];
    return d.norm() > 1e-9 ? d.normalized() : Vec2d(0, 1);
}

static Vec2d exitTangent(const Lane *l) {
    if (!l || l->geometry.points.size() < 2)return Vec2d(0, 1);
    auto &p = l->geometry.points;
    Vec2d d = p[1] - p[0];
    return d.norm() > 1e-9 ? d.normalized() : Vec2d(0, 1);
}

// ─── Group axis helpers ───────────────────────────────────────────────────────
static Vec2d meanPt(const std::vector<ConnId> &cids,
                    const std::vector<Connectivity> &conns,
                    const std::vector<Lane> &lanes, bool use_entry) {
    Vec2d s(0, 0);
    int n = 0;
    for (auto &cid:cids) {
        auto *c = findConn(conns, cid);
        if (!c) continue;
        const Lane *l = use_entry ? findLane(lanes, c->entry_lane_id) : findLane(lanes, c->exit_lane_id);
        s += use_entry ? entryEndpoint(l) : exitEndpoint(l);
        ++n;
    }
    return n > 0 ? (s / n) : Vec2d(0, 0);
}

static Vec2d armLeftNormal(const std::vector<ConnId> &cids,
                           const std::vector<Connectivity> &conns,
                           const std::vector<Lane> &lanes, bool use_entry) {
    Vec2d mt(0, 0);
    int n = 0;
    for (auto &cid:cids) {
        auto *c = findConn(conns, cid);
        if (!c) continue;
        const Lane *l = use_entry ? findLane(lanes, c->entry_lane_id) : findLane(lanes, c->exit_lane_id);
        mt += use_entry ? entryTangent(l) : exitTangent(l);
        ++n;
    }
    if (n == 0)return Vec2d(0, 1);
    mt /= n;
    if (mt.norm() < 1e-9)mt = Vec2d(0, 1);
    mt.normalize();
    return Vec2d(-mt[1], mt[0]);
}

// ─── Pair-specific ref_perp ───────────────────────────────────────────────────
static Vec2d pairRefPerp(const Connectivity &ca, const Connectivity &cb,
                         const std::vector<Lane> &lanes) {
    Vec2d me = 0.5 * (entryEndpoint(findLane(lanes, ca.entry_lane_id))
                      + entryEndpoint(findLane(lanes, cb.entry_lane_id)));
    Vec2d mx = 0.5 * (exitEndpoint(findLane(lanes, ca.exit_lane_id))
                      + exitEndpoint(findLane(lanes, cb.exit_lane_id)));
    Vec2d al = mx - me;
    double n = al.norm();
    if (n < 1e-9)return Vec2d(0, 1);
    al.normalize();
    return Vec2d(-al[1], al[0]);
}

// ─── Sort key computations ────────────────────────────────────────────────────
// Entry cluster: angle of exit_pt relative to mean_entry_pt (DESCENDING = CCW first = LEFT)
static double entryAngle(const Connectivity &c, const Vec2d &me, const std::vector<Lane> &lanes) {
    Vec2d diff = exitEndpoint(findLane(lanes, c.exit_lane_id)) - me;
    return std::atan2(diff[1], diff[0]);
}

// Exit cluster: projection of entry_pt onto exit arm left-normal (DESCENDING = LEFT first)
static double exitEntryLat(const Connectivity &c, const Vec2d &ealn, const std::vector<Lane> &lanes) {
    return entryEndpoint(findLane(lanes, c.entry_lane_id)).dot(ealn);
}

// ─── Pair management ─────────────────────────────────────────────────────────
bool ClusterOrderSolver::hasPair(
    const std::vector<CurvePair>& pairs, const ConnId& a, const ConnId& b) {
    for (auto& p : pairs)
        if ((p.id_a == a && p.id_b == b) || (p.id_a == b && p.id_b == a))
            return true;
    return false;
}

void ClusterOrderSolver::addPairsFromSortedCluster(
        const std::vector<ConnId> &sorted_ids,
        const std::vector<Connectivity> &conns,
        const std::vector<Lane> &lanes) {
    for (int i = 0; i < (int) sorted_ids.size(); ++i) {
        for (int j = i + 1; j < (int) sorted_ids.size(); ++j) {
            if (hasPair(pairs_, sorted_ids[i], sorted_ids[j]))continue;
            auto *ca = findConn(conns, sorted_ids[i]);
            auto *cb = findConn(conns, sorted_ids[j]);
            if (!ca || !cb)continue;
            CurvePair p;
            p.id_a = sorted_ids[i];
            p.id_b = sorted_ids[j];
            p.expected_side = +1;           // sorted[i] is LEFT of sorted[j]
            p.ref_perp = pairRefPerp(*ca, *cb, lanes);

            // ── Shared exit-lane handling ──────────────────────────────────────
            if (!ca->exit_lane_id.empty() &&
                ca->exit_lane_id == cb->exit_lane_id) {
                if (ca->enterGroupId == cb->enterGroupId) {
                    // Same entry group, same exit lane: ordering is still
                    // meaningful throughout most of the path.  Use a LARGER
                    // endpoint skip zone (25%) rather than marking StructuralCross.
                    p.exempt = CrossExemption::None;
                    p.shared_endpoint = true;
                } else {
                    // Different entry groups, same exit lane: paths come from
                    // opposite sides and MUST cross → StructuralCross.
                    p.exempt = CrossExemption::StructuralCross;
                }
                pairs_.push_back(p);
                continue;
            }

            // ── Shared entry-lane handling ─────────────────────────────────────
            // Two curves sharing the same ENTRY lane start from the same (x,y)
            // with the same tangent direction.  The cluster constraint near the
            // shared start would generate a near-zero gradient; use a wider skip
            // zone so the ordering is only enforced in the diverged portion.
            if (!ca->entry_lane_id.empty() &&
                ca->entry_lane_id == cb->entry_lane_id) {
                p.exempt = CrossExemption::None;
                p.shared_endpoint = true;
            } else {
                p.exempt = CrossExemption::None;
            }
            pairs_.push_back(p);
        }
    }
}

// ─── Topological inversion detector ──────────────────────────────────────────
//
// ── Original logic (kept for non-U-turn same-arm pairs) ──────────────────────
//   For entry-cluster pair (A,B):
//     rank_diff = rank_A - rank_B  (sort: DESCENDING angle → rank 0 = most LEFT)
//     lat_diff  = exit_lat_A - exit_lat_B  (in entry-arm perp frame)
//     Consistent: rank_diff * lat_diff < 0  → no structural cross needed
//     Inverted:   rank_diff * lat_diff > 0  → StructuralCross
//
//   For exit-cluster pair (A,B):
//     rank_diff = rank_A - rank_B  (sort: DESCENDING entry_lat → rank 0 = most LEFT)
//     lat_diff  = exit_lat_in_exit_group_A - exit_lat_in_exit_group_B
//                 (exit endpoint projected onto exit-arm left-normal, relative to
//                  mean exit point — INDEPENDENT of the sort key, unlike the old
//                  entry_lat which is identical to the sort key and always gives a
//                  negative product → StructuralCross never fired).
//     Consistent: rank_diff * lat_diff < 0  → no structural cross needed
//     Inverted:   rank_diff * lat_diff > 0  → StructuralCross
void ClusterOrderSolver::detectTopologicalInversions(
        const std::vector<Connectivity> &conns) {
    constexpr double LAT_EPS = 0.3; // m — ignore near-equal laterals
    constexpr double OPP_SIDE_THRESH = 3.0;  // m — minimum lat magnitude for opposite-side check

    for (auto &pair:pairs_) {
        if (pair.exempt != CrossExemption::None)continue;
        auto *ca = findConn(conns, pair.id_a);
        auto *cb = findConn(conns, pair.id_b);
        if (!ca || !cb)continue;

        bool same_entry = (ca->enterGroupId == cb->enterGroupId);
        bool same_exit = (ca->exitGroupId == cb->exitGroupId);

        // ── U-turn vs non-U-turn always produces a structural cross ──
        // A U-turn arc (t0·t1 < -0.5) from the same entry/exit group as a
        // non-U-turn curve MUST pass over that curve — mandatory structural cross.
        // Original code used enterGroupId == exitGroupId which fails when the two
        // group IDs differ even though the curve physically reverses direction.
        {
            bool a_ut = is_uturn_.count(pair.id_a) && is_uturn_.at(pair.id_a);
            bool b_ut = is_uturn_.count(pair.id_b) && is_uturn_.at(pair.id_b);
            if (a_ut != b_ut) {
                pair.exempt = CrossExemption::StructuralCross;
                continue;
            }
        }

        if (same_entry) {
            auto ira = entry_cluster_rank_.find(pair.id_a);
            auto irb = entry_cluster_rank_.find(pair.id_b);
            if (ira == entry_cluster_rank_.end() || irb == entry_cluster_rank_.end()) continue;
            double rank_diff = (double)(ira->second - irb->second);

            // ── U-turn apex inversion ─────────────────────────────────────────
            // For two U-turns from the same entry group, the apex height is
            // proportional to Δlat = |entry_lat − exit_lat| (how much the curve
            // must travel laterally to reach its exit).  When the MORE-RIGHT entry
            // curve (higher rank) also has a LARGER Δlat, its apex swings higher,
            // INVERTING the entry ordering at the apex → unavoidable structural
            // cross.  Detect via rank_diff × dlat_diff > 0.
            {
                bool a_ut = is_uturn_.count(pair.id_a) && is_uturn_.at(pair.id_a);
                bool b_ut = is_uturn_.count(pair.id_b) && is_uturn_.at(pair.id_b);
                if (a_ut && b_ut) {
                    auto iea = entry_lat_in_entry_ref_.find(pair.id_a);
                    auto ieb = entry_lat_in_entry_ref_.find(pair.id_b);
                    auto ixa = exit_lat_in_entry_ref_.find(pair.id_a);
                    auto ixb = exit_lat_in_entry_ref_.find(pair.id_b);
                    if (iea != entry_lat_in_entry_ref_.end() &&
                        ieb != entry_lat_in_entry_ref_.end() &&
                        ixa != exit_lat_in_entry_ref_.end() &&
                        ixb != exit_lat_in_entry_ref_.end()) {
                        double dlat_a = std::abs(ixa->second - iea->second);
                        double dlat_b = std::abs(ixb->second - ieb->second);
                        double dlat_diff = dlat_a - dlat_b;
                        if (std::abs(dlat_diff) > LAT_EPS && rank_diff * dlat_diff > 0.0) {
                            pair.exempt = CrossExemption::StructuralCross;
                            continue;
                        }
                    }
                }
            }

            // ── opposite-side exit = always structural cross ──────────
            {
                auto ixa = exit_lat_sign_.find(pair.id_a);
                auto ixb = exit_lat_sign_.find(pair.id_b);
                if (ixa != exit_lat_sign_.end() && ixb != exit_lat_sign_.end()) {
                    double la = ixa->second, lb = ixb->second;
                    if (la * lb < 0.0 &&
                        std::abs(la) > OPP_SIDE_THRESH &&
                        std::abs(lb) > OPP_SIDE_THRESH) {
                        pair.exempt = CrossExemption::StructuralCross;
                        continue;
                    }
                }
            }

            // ── entry-exit lateral inversion ───────────────────────────
            // Catches cases where the within-arm lane ORDER at entry is INVERTED
            // relative to lane order at exit (e.g. outer right turn from inner lane
            // must cross inner right turn from outer lane; straight east vs right
            // turn south from same west arm).
            // Criterion: (entry_lat_A − entry_lat_B) × (exit_lat_A − exit_lat_B) < 0
            {
                auto iea = entry_lat_in_entry_ref_.find(pair.id_a);
                auto ieb = entry_lat_in_entry_ref_.find(pair.id_b);
                bool share_exit = (!ca->exitGroupId.empty() &&
                                   ca->exitGroupId == cb->exitGroupId);
                if (iea != entry_lat_in_entry_ref_.end() &&
                    ieb != entry_lat_in_entry_ref_.end()) {
                    double entry_diff = iea->second - ieb->second;
                    if (std::abs(entry_diff) > LAT_EPS) {
                        double exit_diff = 0.0;
                        bool got_exit = false;
                        if (share_exit) {
                            auto ixa = exit_lat_in_exit_group_ref_.find(pair.id_a);
                            auto ixb = exit_lat_in_exit_group_ref_.find(pair.id_b);
                            if (ixa != exit_lat_in_exit_group_ref_.end() &&
                                ixb != exit_lat_in_exit_group_ref_.end()) {
                                exit_diff = ixa->second - ixb->second;
                                got_exit = true;
                            }
                        } else {
                            auto ixa = exit_lat_sign_.find(pair.id_a);
                            auto ixb = exit_lat_sign_.find(pair.id_b);
                            if (ixa != exit_lat_sign_.end() &&
                                ixb != exit_lat_sign_.end()) {
                                exit_diff = ixa->second - ixb->second;
                                got_exit = true;
                            }
                        }
                        if (got_exit && std::abs(exit_diff) > LAT_EPS) {
                            if (entry_diff * exit_diff < 0.0) {
                                pair.exempt = CrossExemption::StructuralCross;
                                continue;
                            }
                        }
                    }
                }
            }

            // ── Bug4: same-exit-group pairs — use exit-group lateral ─────────────
            bool share_exit = (!ca->exitGroupId.empty() && ca->exitGroupId == cb->exitGroupId);
            double lat_diff;
            if (share_exit) {
                auto ixa = exit_lat_in_exit_group_ref_.find(pair.id_a);
                auto ixb = exit_lat_in_exit_group_ref_.find(pair.id_b);
                if (ixa == exit_lat_in_exit_group_ref_.end() ||
                    ixb == exit_lat_in_exit_group_ref_.end()) continue;
                lat_diff = ixa->second - ixb->second;
            } else {
                auto ixa = exit_lat_in_entry_ref_.find(pair.id_a);
                auto ixb = exit_lat_in_entry_ref_.find(pair.id_b);
                if (ixa == exit_lat_in_entry_ref_.end() ||
                    ixb == exit_lat_in_entry_ref_.end()) continue;
                lat_diff = ixa->second - ixb->second;
            }
            if (std::abs(lat_diff) < LAT_EPS) continue;
            if (rank_diff * lat_diff > 0)
                pair.exempt = CrossExemption::StructuralCross;

        } else if (same_exit) {
            // use exit-group lateral instead of entry_lat
            // entry_lat_in_exit_ref_ == sort key → product always < 0 → never fires.
            auto ira = exit_cluster_rank_.find(pair.id_a);
            auto irb = exit_cluster_rank_.find(pair.id_b);
            auto ixa = exit_lat_in_exit_group_ref_.find(pair.id_a);
            auto ixb = exit_lat_in_exit_group_ref_.find(pair.id_b);
            if (ira == exit_cluster_rank_.end() || irb == exit_cluster_rank_.end()) continue;
            if (ixa == exit_lat_in_exit_group_ref_.end() ||
                ixb == exit_lat_in_exit_group_ref_.end()) continue;
            double rank_diff = (double)(ira->second - irb->second);
            double lat_diff  = ixa->second - ixb->second;
            if (std::abs(lat_diff) < LAT_EPS) continue;
            if (rank_diff * lat_diff > 0)
                pair.exempt = CrossExemption::StructuralCross;
        }
    }
}

// ─── Main build ───────────────────────────────────────────────────────────────
void ClusterOrderSolver::build(
    const std::vector<Connectivity>& conns,
    const std::vector<Lane>& lanes, const std::vector<LaneGroup>& laneGroups) {
    entry_group_order_.clear();
    exit_group_order_.clear();
    pairs_.clear();
    is_uturn_.clear();
    exit_lat_sign_.clear();
    entry_lat_in_entry_ref_.clear();
    exit_lat_in_entry_ref_.clear();
    entry_lat_in_exit_ref_.clear();
    exit_lat_in_exit_group_ref_.clear();
    entry_cluster_rank_.clear();
    exit_cluster_rank_.clear();

    // Step 1: group by enterGroupId / exitGroupId
    for (auto& c : conns) {
        if (!c.enterGroupId.empty())entry_group_order_[c.enterGroupId].push_back(c.id);
        if (!c.exitGroupId.empty())exit_group_order_[c.exitGroupId].push_back(c.id);
    }

    // Step 2: pre-compute exit-lane lateral positions within each exit group
    // Used as stable tie-break when entry-cluster angles are nearly equal (same direction turns)
    std::unordered_map<LaneId, double> exit_lane_lat_in_group;
    for (auto &kv:exit_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() < 2)continue;
        Vec2d perp = armLeftNormal(cids, conns, lanes, false);
        Vec2d ref = meanPt(cids, conns, lanes, false);
        for (auto &cid:cids) {
            auto *c = findConn(conns, cid);
            if (!c)continue;
            const Lane *xl = findLane(lanes, c->exit_lane_id);
            exit_lane_lat_in_group[c->exit_lane_id] = (exitEndpoint(xl) - ref).dot(perp);
        }
    }

    // Step 3: pre-compute entry-lane lateral positions within each entry group
    std::unordered_map<LaneId, double> entry_lane_lat_in_group;
    for (auto &kv:entry_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() < 2)continue;
        Vec2d perp = armLeftNormal(cids, conns, lanes, true);
        Vec2d ref = meanPt(cids, conns, lanes, true);
        for (auto &cid:cids) {
            auto *c = findConn(conns, cid);
            if (!c)continue;
            const Lane *el = findLane(lanes, c->entry_lane_id);
            entry_lane_lat_in_group[c->entry_lane_id] = (entryEndpoint(el) - ref).dot(perp);
        }
    }

    // Step 4: sort entry clusters
    //   Primary:   DESCENDING angle of (exit_pt - mean_entry_pt)  — global L→R ordering
    //   Tie-break: DESCENDING exit_lane lateral position in exit group — stable same-dir ordering
    //   Final:     DESCENDING entry_lane lateral position in entry group
    for (auto &kv:entry_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() < 2)continue;
        Vec2d me = meanPt(cids, conns, lanes, true);
        std::stable_sort(cids.begin(), cids.end(), [&](const ConnId &a, const ConnId &b) {
            auto *ca = findConn(conns, a);
            auto *cb = findConn(conns, b);
            if (!ca || !cb)return false;
            double ang_a = entryAngle(*ca, me, lanes);
            double ang_b = entryAngle(*cb, me, lanes);
            if (std::abs(ang_a - ang_b) > 0.035)return ang_a > ang_b; // >2° angle decides
            // Same-direction multi-lane: use exit group lateral ordering
            double xl_a = exit_lane_lat_in_group.count(ca->exit_lane_id) ?
                          exit_lane_lat_in_group.at(ca->exit_lane_id) : 0;
            double xl_b = exit_lane_lat_in_group.count(cb->exit_lane_id) ?
                          exit_lane_lat_in_group.at(cb->exit_lane_id) : 0;
            if (std::abs(xl_a - xl_b) > 0.1)return xl_a > xl_b;
            // Entry-lane lateral tie-break
            double el_a = entry_lane_lat_in_group.count(ca->entry_lane_id) ?
                          entry_lane_lat_in_group.at(ca->entry_lane_id) : 0;
            double el_b = entry_lane_lat_in_group.count(cb->entry_lane_id) ?
                          entry_lane_lat_in_group.at(cb->entry_lane_id) : 0;
            return el_a > el_b;
        });
    }

    // Step 5: sort exit clusters
    //   Primary:   DESCENDING entry_lat on exit arm left-normal
    //   Tie-break: DESCENDING exit_lane lateral position in exit group
    for (auto &kv:exit_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() < 2)continue;
        Vec2d ealn = armLeftNormal(cids, conns, lanes, false);
        std::stable_sort(cids.begin(), cids.end(), [&](const ConnId &a, const ConnId &b) {
            auto *ca = findConn(conns, a);
            auto *cb = findConn(conns, b);
            if (!ca || !cb)return false;
            double la = exitEntryLat(*ca, ealn, lanes);
            double lb = exitEntryLat(*cb, ealn, lanes);
            if (std::abs(la - lb) > 0.5)return la > lb;
            double xl_a = exit_lane_lat_in_group.count(ca->exit_lane_id) ?
                          exit_lane_lat_in_group.at(ca->exit_lane_id) : 0;
            double xl_b = exit_lane_lat_in_group.count(cb->exit_lane_id) ?
                          exit_lane_lat_in_group.at(cb->exit_lane_id) : 0;
            return xl_a > xl_b;
        });
    }

    // Step 6: store ranks and lateral positions for topological inversion detection
    // Also pre-compute U-turn flags via actual t0·t1 dot products.
    for (auto &kv:entry_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        Vec2d me = meanPt(cids, conns, lanes, true);
        Vec2d perp = armLeftNormal(cids, conns, lanes, true);
        for (int i = 0; i < (int) cids.size(); ++i) {
            entry_cluster_rank_[cids[i]] = i;
            auto *c = findConn(conns, cids[i]);
            if (!c)continue;
            const Lane *el = findLane(lanes, c->entry_lane_id);
            const Lane *xl = findLane(lanes, c->exit_lane_id);
            // exit_lat in entry arm reference (higher = more LEFT)
            double elat = (exitEndpoint(xl) - me).dot(perp);
            exit_lat_in_entry_ref_[cids[i]] = elat;
            exit_lat_sign_[cids[i]] = elat; // raw value for sign comparison

            // entry endpoint lateral within its own entry-group frame.
            const Lane *entl = findLane(lanes, c->entry_lane_id);
            entry_lat_in_entry_ref_[cids[i]] = (entryEndpoint(entl) - me).dot(perp);
            // detect U-turns by t0·t1 dot product using actual lane tangents
            // enterGroupId == exitGroupId check was wrong: physical U-turns can have
            // different group IDs (one entry group, one exit group, same physical arm).
            Vec2d t0 = entryTangent(el);
            Vec2d t1 = exitTangent(xl);
            is_uturn_[cids[i]] = (t0.dot(t1) < -0.5);
        }
    }

    for (auto &kv:exit_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        Vec2d ealn = armLeftNormal(cids, conns, lanes, false);
        // compute exit group's own reference for exit lateral positions.
        // entry_lat_in_exit_ref_ is the SAME metric used for sorting → rank_diff and
        // lat_diff always have opposite signs → StructuralCross never fires for same_exit.
        // exit_lat_in_exit_group_ref_ uses the EXIT endpoint in the exit-arm frame,
        // which is independent of the sort key and can reveal true geometric inversions.
        Vec2d exit_ref_pt = meanPt(cids, conns, lanes, false);
        for (int i = 0; i < (int) cids.size(); ++i) {
            exit_cluster_rank_[cids[i]] = i;
            auto *c = findConn(conns, cids[i]);
            if (!c)continue;
            const Lane *el = findLane(lanes, c->entry_lane_id);
            const Lane *xl = findLane(lanes, c->exit_lane_id);
            // entry_lat in exit arm reference (used for existing same_exit check — kept for compat.)
            entry_lat_in_exit_ref_[cids[i]] = entryEndpoint(el).dot(ealn);
            // exit_lat in exit group's own frame: (exit_endpoint − mean_exit) · exit_arm_left_normal
            exit_lat_in_exit_group_ref_[cids[i]] = (exitEndpoint(xl) - exit_ref_pt).dot(ealn);
            // Also populate is_uturn_ for conns not appearing in any entry group
            if (!is_uturn_.count(cids[i])) {
                Vec2d t0 = entryTangent(el);
                Vec2d t1 = exitTangent(xl);
                is_uturn_[cids[i]] = (t0.dot(t1) < -0.5);
            }
        }
    }

    // Step 7: build constraint pairs from entry clusters
    for (auto &kv:entry_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() > 1)addPairsFromSortedCluster(cids, conns, lanes);
    }

    // Step 8: build constraint pairs from exit clusters (new pairs only)
    for (auto &kv:exit_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() > 1) {
            addPairsFromSortedCluster(cids, conns, lanes);
        }
    }

    // Step 9: detect and mark true topological inversions as StructuralCross
    detectTopologicalInversions(conns);
}

// ─── Lookup helpers ───────────────────────────────────────────────────────────
bool ClusterOrderSolver::pairExists(const ConnId& a, const ConnId& b) const {
    return hasPair(pairs_,a,b);
}

CrossExemption ClusterOrderSolver::exemptionOf(const ConnId& a, const ConnId& b) const {
    for (auto& p : pairs_)
        if ((p.id_a == a && p.id_b == b) || (p.id_a == b && p.id_b == a))
            return p.exempt;
    return CrossExemption::None;
}

int ClusterOrderSolver::expectedSideOf(const ConnId &a, const ConnId &b) const {
    for (auto &p:pairs_) {
        if (p.id_a == a && p.id_b == b)return p.expected_side;
        if (p.id_a == b && p.id_b == a)return -p.expected_side;
    }
    return 0;
}

Vec2d ClusterOrderSolver::refPerpOf(const ConnId &a, const ConnId &b) const {
    for (auto &p:pairs_)
        if ((p.id_a == a && p.id_b == b) || (p.id_a == b && p.id_b == a))return p.ref_perp;
    return Vec2d(0, 0);
}

void ClusterOrderSolver::markObstacleExempt(CurvePair &pair, const Vec2d &pt,
                                            const SDFField &sdf, double r) {
    auto kv = sdf.queryWithGrad(pt);
    auto d = kv.first;
    auto _g = kv.second;
    if (d < r) {
        pair.exempt = CrossExemption::ObstacleCross;
        pair.exempt_zone_radius = r;
    }
}

void ClusterOrderSolver::checkAndMarkA2(
    const std::unordered_map<ConnId, BezierCurve>& curves, const SDFField& sdf, double r) {
    for (auto& pair : pairs_) {
        if (pair.exempt != CrossExemption::None)
            continue;
        auto ia = curves.find(pair.id_a), ib = curves.find(pair.id_b);
        if (ia == curves.end() || ib == curves.end())
            continue;
        for (auto& pt : curveCrossings(ia->second, ib->second))
            markObstacleExempt(pair, pt, sdf, r);
    }
}

}