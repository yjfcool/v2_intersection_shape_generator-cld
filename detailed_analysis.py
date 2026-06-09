#!/usr/bin/env python3
"""
Detailed analysis of the specific issues mentioned in the intersection shape generator:
- U-turn formations conn20|conn22 and conn19|conn21 not forming U-turn shapes
- Non-endpoint intersections in same-cluster connections
- Analysis of cluster sorting problems
"""

import json
import math
from typing import Dict, List, Tuple, Optional


class DetailedConnectionAnalyzer:
    def __init__(self, json_file_path: str):
        """Initialize with the intersection data."""
        with open(json_file_path, 'r', encoding='utf-8') as f:
            self.data = json.load(f)

        self.connections = {c['id']: c for c in self.data['connectivities']}
        self.lanes = {l['id']: l for l in self.data['lanes']}

    def get_lane_endpoint(self, lane_id: str, entry: bool = True) -> Optional[Tuple[float, float]]:
        """Get endpoint of a lane."""
        if lane_id not in self.lanes:
            return None

        lane = self.lanes[lane_id]
        if not lane.get('geometry', {}).get('points'):
            return None

        points = lane['geometry']['points']
        if entry:  # entry is the last point
            return (points[-1]['x'], points[-1]['y'])
        else:  # exit is the first point
            return (points[0]['x'], points[0]['y'])

    def get_lane_tangent(self, lane_id: str, entry: bool = True) -> Optional[Tuple[float, float]]:
        """Get tangent direction of a lane at the specified end."""
        if lane_id not in self.lanes:
            return None

        lane = self.lanes[lane_id]
        if not lane.get('geometry', {}).get('points'):
            return None

        points = lane['geometry']['points']
        if len(points) < 2:
            return None

        if entry:  # entry tangent uses last two points
            dx = points[-1]['x'] - points[-2]['x']
            dy = points[-1]['y'] - points[-2]['y']
        else:  # exit tangent uses first two points
            dx = points[1]['x'] - points[0]['x']
            dy = points[1]['y'] - points[0]['y']

        # Normalize
        mag = math.sqrt(dx*dx + dy*dy)
        if mag < 1e-9:
            return (0.0, 1.0)  # default direction

        return (dx/mag, dy/mag)

    def calculate_angle_between_vectors(self, v1: Tuple[float, float], v2: Tuple[float, float]) -> float:
        """Calculate angle between two vectors in radians."""
        dot = v1[0]*v2[0] + v1[1]*v2[1]
        det = v1[0]*v2[1] - v1[1]*v2[0]
        return math.atan2(det, dot)

    def analyze_connection_geometry(self, conn_id: str):
        """Analyze geometry of a single connection."""
        if conn_id not in self.connections:
            return None

        conn = self.connections[conn_id]

        # Get endpoints and tangents
        entry_pt = self.get_lane_endpoint(conn['entry_lane_id'], entry=True)
        exit_pt = self.get_lane_endpoint(conn['exit_lane_id'], entry=False)
        entry_tan = self.get_lane_tangent(conn['entry_lane_id'], entry=True)
        exit_tan = self.get_lane_tangent(conn['exit_lane_id'], entry=False)

        if not all([entry_pt, exit_pt, entry_tan, exit_tan]):
            return None

        # Calculate entry-exit angle for U-turn check
        entry_exit_angle = self.calculate_angle_between_vectors(entry_tan, exit_tan)

        # Calculate geometric turn type
        path_vec = (exit_pt[0] - entry_pt[0], exit_pt[1] - entry_pt[1])
        path_mag = math.sqrt(path_vec[0]**2 + path_vec[1]**2)

        if path_mag > 1e-9:
            path_dir = (path_vec[0]/path_mag, path_vec[1]/path_mag)
        else:
            path_dir = (1, 0)

        # Cross product to determine turn direction (positive = left turn, negative = right)
        cross = entry_tan[0]*path_dir[1] - entry_tan[1]*path_dir[0]

        if abs(cross) < 0.25:  # Within ~15 degrees of straight
            turn_type = "straight"
        elif cross > 0:
            turn_type = "left"
        else:
            turn_type = "right"

        is_uturn = abs(entry_exit_angle) > math.pi * 0.85  # More than 153 degrees (almost 180)

        return {
            'conn_id': conn_id,
            'entry_point': entry_pt,
            'exit_point': exit_pt,
            'entry_tangent': entry_tan,
            'exit_tangent': exit_tan,
            'entry_exit_angle_deg': math.degrees(entry_exit_angle),
            'geometric_turn_type': turn_type,
            'is_uturn_geometric': is_uturn,
            'declared_turn_type': conn.get('turn_type'),
            'enter_group_id': conn.get('enterGroupId'),
            'exit_group_id': conn.get('exitGroupId')
        }

    def analyze_cluster_sorting_logic(self):
        """Analyze the cluster sorting logic based on the actual code implementation."""
        print("=" * 80)
        print("DETAILED CLUSTER SORTING ANALYSIS")
        print("=" * 80)
        print("\nAccording to the code in src/constraints/cluster_order.cpp:")
        print("1. Entry clusters: Sorted by DESCENDING angle of (exit_pt - mean_entry_pt)")
        print("2. Exit clusters: Sorted by DESCENDING projection of entry_pt onto exit arm left-normal")
        print("3. For same-cluster pairs: Expected side relationships maintained")

        # Group by entry/exit groups
        entry_groups = {}
        exit_groups = {}

        for conn_id, conn in self.connections.items():
            enter_gid = conn.get('enterGroupId', '')
            exit_gid = conn.get('exitGroupId', '')

            if enter_gid:
                if enter_gid not in entry_groups:
                    entry_groups[enter_gid] = []
                entry_groups[enter_gid].append(conn_id)

            if exit_gid:
                if exit_gid not in exit_groups:
                    exit_groups[exit_gid] = []
                exit_groups[exit_gid].append(conn_id)

        print(f"\nENTRY CLUSTERS FOUND ({len(entry_groups)} total):")
        for gid, conn_list in entry_groups.items():
            print(f"  Group {gid}: {sorted(conn_list)}")
            if any(int(x) in [20, 22, 19, 21] for x in conn_list):
                print(f"    → Contains problematic U-turn candidates: {[x for x in conn_list if x in ['20', '22', '19', '21']]}")

        print(f"\nEXIT CLUSTERS FOUND ({len(exit_groups)} total):")
        for gid, conn_list in exit_groups.items():
            print(f"  Group {gid}: {sorted(conn_list)}")
            if any(int(x) in [10, 12, 30, 32, 18] for x in conn_list):
                print(f"    → Contains problematic same-cluster candidates: {[x for x in conn_list if x in ['10', '12', '30', '32', '18']]}")


    def analyze_problematic_connections(self):
        """Analyze the specific problematic connections mentioned in the issue."""
        print("\n" + "=" * 80)
        print("ANALYSIS OF SPECIFIC PROBLEMATIC CONNECTIONS")
        print("=" * 80)

        # U-turn pairs mentioned: conn20|conn22 and conn19|conn21
        uturn_pairs = [('20', '22'), ('19', '21')]

        print("\n1. U-TURN FORMATION ISSUES:")
        for conn1_id, conn2_id in uturn_pairs:
            print(f"\n   U-turn pair: {conn1_id} | {conn2_id}")

            conn1_geo = self.analyze_connection_geometry(conn1_id)
            conn2_geo = self.analyze_connection_geometry(conn2_id)

            if conn1_geo and conn2_geo:
                print(f"     Conn {conn1_id}: Angle={conn1_geo['entry_exit_angle_deg']:.1f}°, "
                      f"Geometric type={conn1_geo['geometric_turn_type']}, "
                      f"Is U-turn: {conn1_geo['is_uturn_geometric']}")

                print(f"     Conn {conn2_id}: Angle={conn2_geo['entry_exit_angle_deg']:.1f}°, "
                      f"Geometric type={conn2_geo['geometric_turn_type']}, "
                      f"Is U-turn: {conn2_geo['is_uturn_geometric']}")

                # Check if they belong to same entry group (which they should based on our cluster analysis)
                same_entry = conn1_geo['enter_group_id'] == conn2_geo['enter_group_id']
                same_exit = conn1_geo['exit_group_id'] == conn2_geo['exit_group_id']

                print(f"     Same entry group: {same_entry} ({conn1_geo['enter_group_id']})")
                print(f"     Same exit group: {same_exit} ({conn1_geo['exit_group_id']})")

                # Root cause analysis
                if (conn1_geo['is_uturn_geometric'] or conn2_geo['is_uturn_geometric']) and \
                   not (conn1_geo['is_uturn_geometric'] and conn2_geo['is_uturn_geometric']):
                    print(f"     ⚠️  ISSUE: Mixed U-turn/non-U-turn in same cluster - this can cause non-U-turn shapes")

        # Same-cluster pairs mentioned: conn10|conn12, conn30|conn12, conn32|conn18
        # Since we saw earlier that conn10 might not exist, let's check available pairs
        print(f"\n2. SAME-CLUSTER NON-ENDPOINT INTERSECTION ISSUES:")

        # First, identify connections that should be in same clusters
        same_entry_pairs = []
        same_exit_pairs = []

        # Go through all connection pairs and identify same-cluster pairs
        for conn1_id in self.connections:
            for conn2_id in self.connections:
                if conn1_id >= conn2_id:  # avoid duplicates
                    continue

                conn1 = self.connections[conn1_id]
                conn2 = self.connections[conn2_id]

                same_entry = conn1.get('enterGroupId') == conn2.get('enterGroupId')
                same_exit = conn1.get('exitGroupId') == conn2.get('exitGroupId')

                if same_entry:
                    same_entry_pairs.append((conn1_id, conn2_id))
                if same_exit:
                    same_exit_pairs.append((conn1_id, conn2_id))

        print(f"\n   Same-entry cluster pairs found: {len(same_entry_pairs)}")
        for pair in same_entry_pairs[:10]:  # Show first 10
            print(f"     {pair[0]} | {pair[1]}")

        print(f"\n   Same-exit cluster pairs found: {len(same_exit_pairs)}")
        for pair in same_exit_pairs[:10]:  # Show first 10
            print(f"     {pair[0]} | {pair[1]}")


    def analyze_root_causes(self):
        """Analyze the root causes of the issues based on code review."""
        print("\n" + "=" * 80)
        print("ROOT CAUSE ANALYSIS")
        print("=" * 80)

        print("\nA. U-TURN FORMATION PROBLEMS:")
        print("   1. The U-turn formation is handled by buildTwoSegmentUTurn() in hermite_init.cpp")
        print("   2. This function expects anti-parallel tangents (dot product < -0.5) to detect U-turns")
        print("   3. Issues arise when:")
        print("      - The entry and exit tangents aren't properly anti-parallel")
        print("      - Obstacle interference forces the curve to deviate from U-turn shape")
        print("      - Lateral offset calculation doesn't account for proper U-turn geometry")

        print("\nB. NON-ENDPOINT INTERSECTION PROBLEMS:")
        print("   1. The cluster ordering logic in ClusterOrderSolver::build() sorts connections within clusters")
        print("   2. Same-entry clusters sorted by descending angle of (exit_pt - mean_entry_pt)")
        print("   3. Same-exit clusters sorted by descending projection of entry_pt onto exit arm left-normal")
        print("   4. Problems occur when:")
        print("      - Initial curve placement doesn't respect the intended ordering")
        print("      - The needsLateralOffset() function fails to apply sufficient separation")
        print("      - Curves cross during optimization despite ordering constraints")

        print("\nC. SPECIFIC CODE AREAS TO CHECK:")
        print("   1. In connectivity_generator.cpp:")
        print("      - needsLateralOffset() function for initial curve positioning")
        print("      - buildSiblings() function for constraint application")
        print("      - generateOne() for the curve generation workflow")
        print("   2. In cluster_order.cpp:")
        print("      - detectTopologicalInversions() for proper exemption marking")
        print("      - Sort key computations for accurate clustering")
        print("   3. In hermite_init.cpp:")
        print("      - buildTwoSegmentUTurn() for proper U-turn geometry")
        print("      - buildInitialCurve() for obstacle-aware initial placement")

    def analyze_solutions(self):
        """Suggest potential solutions based on root cause analysis."""
        print("\n" + "=" * 80)
        print("POTENTIAL SOLUTIONS")
        print("=" * 80)

        print("\n1. FOR U-TURN FORMATION ISSUES:")
        print("   • Enhance U-turn detection logic to handle near-anti-parallel cases")
        print("   • Improve apex calculation in buildTwoSegmentUTurn() to ensure proper semicircular shape")
        print("   • Add obstacle-aware U-turn placement to avoid interference")

        print("\n2. FOR NON-ENDPOINT INTERSECTION ISSUES:")
        print("   • Strengthen the needsLateralOffset() function to provide adequate separation")
        print("   • Improve the initial curve generation to respect cluster ordering from the start")
        print("   • Enhance the detectTopologicalInversions() logic to better identify legitimate crossings")

        print("\n3. GENERAL IMPROVEMENTS:")
        print("   • Add more robust collision detection during curve optimization")
        print("   • Implement iterative refinement for cluster ordering violations")
        print("   • Enhance the adaptive refinement process to handle tight spaces better")


def main():
    analyzer = DetailedConnectionAnalyzer('datas/intersection_cross.json')

    # Perform detailed analysis
    analyzer.analyze_cluster_sorting_logic()
    analyzer.analyze_problematic_connections()
    analyzer.analyze_root_causes()
    analyzer.analyze_solutions()


if __name__ == "__main__":
    main()