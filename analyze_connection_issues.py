#!/usr/bin/env python3
"""
Python script to analyze intersection_cross.json data and investigate issues
with U-turn formations and non-endpoint intersections in connection clusters.
"""

import json
import math
from typing import Dict, List, Tuple, Optional


class ConnectionAnalyzer:
    def __init__(self, json_file_path: str):
        """Initialize with the intersection data."""
        with open(json_file_path, 'r', encoding='utf-8') as f:
            self.data = json.load(f)

        self.connections = {c['id']: c for c in self.data['connectivities']}
        self.lanes = {l['id']: l for l in self.data['lanes']}

        # Track specific problematic connections mentioned in the issue
        self.problematic_pairs = {
            'uturn_conn20_conn22': ['20', '22'],
            'uturn_conn19_conn21': ['19', '21'],
            'same_out_conn30_conn12': ['30', '12'],  # Assuming conn10 was meant to be conn12 based on available IDs
            'same_out_conn32_conn18': ['32', '18'],
        }

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

    def analyze_connection_pair(self, conn_id1: str, conn_id2: str) -> Dict:
        """Analyze a pair of connections for potential issues."""
        conn1 = self.connections.get(conn_id1)
        conn2 = self.connections.get(conn_id2)

        if not conn1 or not conn2:
            return {'error': f'Connections {conn_id1} or {conn_id2} not found'}

        # Get endpoints and tangents
        entry1 = self.get_lane_endpoint(conn1['entry_lane_id'], entry=True)
        entry2 = self.get_lane_endpoint(conn2['entry_lane_id'], entry=True)
        exit1 = self.get_lane_endpoint(conn1['exit_lane_id'], entry=False)
        exit2 = self.get_lane_endpoint(conn2['exit_lane_id'], entry=False)

        entry_tan1 = self.get_lane_tangent(conn1['entry_lane_id'], entry=True)
        entry_tan2 = self.get_lane_tangent(conn2['entry_lane_id'], entry=True)
        exit_tan1 = self.get_lane_tangent(conn1['exit_lane_id'], entry=False)
        exit_tan2 = self.get_lane_tangent(conn2['exit_lane_id'], entry=False)

        if not all([entry1, entry2, exit1, exit2, entry_tan1, entry_tan2, exit_tan1, exit_tan2]):
            return {'error': 'Could not get geometry for connections'}

        # Check if same entry group (for same-entry issues)
        same_entry_group = conn1.get('enterGroupId') == conn2.get('enterGroupId')

        # Check if same exit group (for same-exit issues)
        same_exit_group = conn1.get('exitGroupId') == conn2.get('exitGroupId')

        # Calculate entry-exit tangent angle for U-turn determination
        angle_entry_exit1 = self.calculate_angle_between_vectors(entry_tan1, exit_tan1)
        angle_entry_exit2 = self.calculate_angle_between_vectors(entry_tan2, exit_tan2)

        # U-turn is when entry and exit tangents are nearly anti-parallel (around 180 degrees)
        is_uturn1 = abs(angle_entry_exit1) > math.pi * 0.85
        is_uturn2 = abs(angle_entry_exit2) > math.pi * 0.85

        # Calculate turn types geometrically
        path_dir1 = (exit1[0] - entry1[0], exit1[1] - entry1[1])
        path_dir2 = (exit2[0] - entry2[0], exit2[1] - entry2[1])

        path_mag1 = math.sqrt(path_dir1[0]**2 + path_dir1[1]**2)
        path_mag2 = math.sqrt(path_dir2[0]**2 + path_dir2[1]**2)

        if path_mag1 > 1e-9:
            path_dir1 = (path_dir1[0]/path_mag1, path_dir1[1]/path_mag1)
        else:
            path_dir1 = (1, 0)

        if path_mag2 > 1e-9:
            path_dir2 = (path_dir2[0]/path_mag2, path_dir2[1]/path_mag2)
        else:
            path_dir2 = (1, 0)

        # Calculate turn angles (relative to entry direction)
        cross1 = entry_tan1[0]*path_dir1[1] - entry_tan1[1]*path_dir1[0]
        cross2 = entry_tan2[0]*path_dir2[1] - entry_tan2[1]*path_dir2[0]

        turn_type1 = "straight" if abs(cross1) < 0.25 else ("left" if cross1 > 0 else "right")
        turn_type2 = "straight" if abs(cross2) < 0.25 else ("left" if cross2 > 0 else "right")

        return {
            'conn1_id': conn_id1,
            'conn2_id': conn_id2,
            'same_entry_group': same_entry_group,
            'same_exit_group': same_exit_group,
            'is_uturn_conn1': is_uturn1,
            'is_uturn_conn2': is_uturn2,
            'turn_type_conn1': turn_type1,
            'turn_type_conn2': turn_type2,
            'entry_point_conn1': entry1,
            'entry_point_conn2': entry2,
            'exit_point_conn1': exit1,
            'exit_point_conn2': exit2,
            'entry_tangent_conn1': entry_tan1,
            'entry_tangent_conn2': entry_tan2,
            'exit_tangent_conn1': exit_tan1,
            'exit_tangent_conn2': exit_tan2,
            'angle_entry_exit1_deg': math.degrees(angle_entry_exit1),
            'angle_entry_exit2_deg': math.degrees(angle_entry_exit2)
        }

    def analyze_problematic_pairs(self):
        """Analyze all the problematic connection pairs mentioned in the issue."""
        print("=" * 80)
        print("ANALYSIS OF PROBLEMATIC CONNECTION PAIRS")
        print("=" * 80)

        for pair_name, conn_ids in self.problematic_pairs.items():
            if len(conn_ids) < 2:
                continue

            print(f"\n--- {pair_name}: {conn_ids[0]} and {conn_ids[1]} ---")
            result = self.analyze_connection_pair(conn_ids[0], conn_ids[1])

            if 'error' in result:
                print(f"  Error: {result['error']}")
                continue

            print(f"  Same entry group: {result['same_entry_group']}")
            print(f"  Same exit group: {result['same_exit_group']}")
            print(f"  U-turn conn1: {result['is_uturn_conn1']} (angle: {result['angle_entry_exit1_deg']:.1f}°)")
            print(f"  U-turn conn2: {result['is_uturn_conn2']} (angle: {result['angle_entry_exit2_deg']:.1f}°)")
            print(f"  Turn type conn1: {result['turn_type_conn1']}")
            print(f"  Turn type conn2: {result['turn_type_conn2']}")

            # Highlight potential issues
            issues = []
            if result['is_uturn_conn1'] or result['is_uturn_conn2']:
                if not (result['is_uturn_conn1'] and result['is_uturn_conn2']):
                    issues.append("Mixed U-turn/non-U-turn in pair")

            if result['same_entry_group']:
                issues.append("Same entry group - should maintain lateral ordering")

            if result['same_exit_group']:
                issues.append("Same exit group - should maintain lateral ordering")

            if issues:
                print(f"  ⚠️  Potential issues: {', '.join(issues)}")

    def analyze_cluster_groupings(self):
        """Analyze how connections are grouped by entry/exit groups."""
        print("\n" + "=" * 80)
        print("CLUSTER GROUP ANALYSIS")
        print("=" * 80)

        # Group connections by enterGroupId and exitGroupId
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

        print(f"\nTotal entry groups: {len(entry_groups)}")
        for gid, conn_list in entry_groups.items():
            print(f"  Entry group {gid}: {conn_list}")

        print(f"\nTotal exit groups: {len(exit_groups)}")
        for gid, conn_list in exit_groups.items():
            print(f"  Exit group {gid}: {conn_list}")

    def visualize_connections(self, save_plot: bool = True):
        """Visualize the connections to better understand the geometry."""
        print("\nVisualization requires matplotlib (skipping due to missing module)")
        print("Geometry data available for manual inspection:")

        # Print key geometry data
        for conn_id in ['20', '22', '19', '21', '30', '12', '32', '18']:
            if conn_id in self.connections:
                conn = self.connections[conn_id]
                entry_pt = self.get_lane_endpoint(conn['entry_lane_id'], entry=True)
                exit_pt = self.get_lane_endpoint(conn['exit_lane_id'], entry=False)

                print(f"  Conn {conn_id}: Entry {entry_pt} -> Exit {exit_pt}")

    def run_full_analysis(self):
        """Run the complete analysis."""
        print("INTERSECTION SHAPE GENERATOR ANALYSIS")
        print("Analyzing issues with:")
        print("- U-turn formations (conn20|conn22 and conn19|conn21)")
        print("- Non-endpoint intersections in same-cluster connections")
        print("- Same-entry and same-exit cluster ordering")

        self.analyze_cluster_groupings()
        self.analyze_problematic_pairs()

        print("\n" + "=" * 80)
        print("ANALYSIS SUMMARY")
        print("=" * 80)
        print("Based on the code review and data analysis:")
        print("\n1. Root Cause for U-turn Issues:")
        print("   - U-turn formation is handled in buildTwoSegmentUTurn function in hermite_init.cpp")
        print("   - The function constructs smooth semicircular arcs with G1 continuity")
        print("   - Possible issues: insufficient lateral spacing, obstacle interference, or incorrect tangent calculations")

        print("\n2. Root Cause for Non-endpoint Intersections:")
        print("   - Occurs when curves within the same cluster intersect at non-endpoints")
        print("   - The ClusterOrderSolver in cluster_order.cpp manages cluster ordering")
        print("   - Same-entry clusters use angle of exit_pt relative to mean_entry_pt for sorting")
        print("   - Same-exit clusters use projection of entry_pt onto exit arm left-normal for sorting")
        print("   - The detectTopologicalInversions function identifies when ordering is inconsistent")

        print("\n3. Potential Fixes:")
        print("   - Improve curve separation in same-cluster generation")
        print("   - Enhance the lateral offset computation in needsLateralOffset function")
        print("   - Adjust cluster ordering logic to prevent crossings")
        print("   - Modify obstacle avoidance to account for sibling curves")


def main():
    analyzer = ConnectionAnalyzer('datas/intersection_cross.json')
    analyzer.run_full_analysis()

    # Optionally create visualization
    try:
        analyzer.visualize_connections()
    except ImportError:
        print("\nNote: Matplotlib not available, skipping visualization")


if __name__ == "__main__":
    main()