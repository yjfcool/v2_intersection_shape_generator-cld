//
// Created by ubuntu on 26-5-24.
//

#ifndef IODATA_SHAPEFILE_H
#define IODATA_SHAPEFILE_H
#pragma once

#include <vector>
#include <iostream>
#include "shapefile.hpp"
#include "filesystem.hpp"
#include "types.h"

namespace isg {

namespace fs = ghc::filesystem;

inline std::vector<ShapePoint> toShapepoints(const std::vector<Vec2d>& points) {
    std::vector<ShapePoint> shppoints;
    for (auto p : points) {
        shppoints.emplace_back(ShapePoint{p[0], p[1], 0.0});
    }
    return shppoints;
};

static bool save(IntersectionInput& inp, std::string out_dir, std::string prefix) {
    auto writeLanes = [&](const std::string& dir, const std::string& fname,
                          const std::vector<Lane>& lanes) {
        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"LANE_ORDER", 'C', 64, 0},
            {"GROUP_ID", 'C', 64, 0},
            {"GROUP_TYPE", 'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < lanes.size(); ++i) {
            const auto& pl = lanes[i];
            std::vector<std::string> attrs = {pl.id, "-1", pl.groupId, ""};
            std::vector<ShapePoint> points = toShapepoints(pl.geometry.points);
            ShapeRecord record(i, shpType, attrs, points, {0});
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeLaneEdges = [&](const std::string& dir, const std::string& fname,
                              const std::vector<LaneEdge>& edgelines) {
        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"GROUP_ID", 'C', 64, 0},
            {"GROUP_TYPE", 'C', 64, 0},
            {"LEFT_CLINE_ID", 'C', 64, 0},
            {"RIGHT_CLINE_ID", 'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < edgelines.size(); ++i) {
            const auto& pl = edgelines[i];
            std::vector<std::string> attrs = {pl.id, "-1", pl.groupId, "", ""};
            std::vector<ShapePoint> points = toShapepoints(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeStopLines = [&](const std::string& dir, const std::string& fname,
                              const std::vector<StopLine>& stoplines) {
        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"ENTRY_GROUP_ID", 'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < stoplines.size(); ++i) {
            const auto& pl = stoplines[i];
            std::vector<std::string> attrs = {pl.id, pl.associated_group_id};
            std::vector<ShapePoint> points = toShapepoints(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeRoadEdges = [&](const std::string& dir, const std::string& fname,
                              const std::vector<Boundary>& boundaries) {
        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"TYPE", 'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < boundaries.size(); ++i) {
            const auto& pl = boundaries[i];
            std::vector<std::string> attrs = {pl.id, std::to_string((int)pl.type)};
            std::vector<ShapePoint> points = toShapepoints(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeObstacles = [&](const std::string& dir, const std::string& fname,
                              const std::vector<Obstacle>& obstacles) {
        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"TYPE", 'C', 64, 0},
        };
        int shpType = SHP_POLYGONZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < obstacles.size(); ++i) {
            const auto& pg = obstacles[i];
            std::vector<std::string> attrs = {pg.id, "-1"}; //std::to_string(pg.type)
            std::vector<ShapePoint> points = toShapepoints(pg.geometry.outer);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeCrosswalks = [&](const std::string& dir, const std::string& fname,
                               const std::vector<Crosswalk>& crosswalks) {
        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"TYPE", 'C', 64, 0},
        };
        int shpType = SHP_POLYGONZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < crosswalks.size(); ++i) {
            const auto& pg = crosswalks[i];
            std::vector<std::string> attrs = {pg.id, "-1"}; //std::to_string(pg.type)
            std::vector<ShapePoint> points = toShapepoints(pg.geometry.outer);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeAreas = [&](const std::string& dir, const std::string& fname,
                          const IntersectionArea& inter_area) {
        std::vector<Polygon2d> areas = {inter_area.geometry};

        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"TYPE", 'C', 64, 0},
        };
        int shpType = SHP_POLYGONZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < areas.size(); ++i) {
            const auto& pg = areas[i];
            std::vector<std::string> attrs = {std::to_string(i), "-1"}; //std::to_string(pg.type)
            std::vector<ShapePoint> points = toShapepoints(pg.outer);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    writeLanes(out_dir, prefix + "_lanes", inp.lanes);
    writeLaneEdges(out_dir, prefix + "_laneedges", inp.lane_edges);
    writeStopLines(out_dir, prefix + "_stoplines", inp.stop_lines);
    writeRoadEdges(out_dir, prefix + "_roadedges", inp.boundaries);
    writeObstacles(out_dir, prefix + "_obstacles", inp.obstacles);
    writeCrosswalks(out_dir, prefix + "_crosswalks", inp.crosswalks);
    writeAreas(out_dir, prefix + "_areas", inp.area);
    return true;
}

static bool save(IntersectionOutput& out, std::string out_dir, std::string prefix) {
    auto writeLanes = [&](const std::string& dir, const std::string& fname,
                          const std::vector<ConnectivityCurve>& lanes) {
        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"TURN_TYPE", 'C', 64, 0},
            {"FLANE", 'C', 64, 0},
            {"TLANE", 'C', 64, 0},
            {"LANE_TYPE", 'C', 64, 0},
            {"FIXED_SHAPE", 'C', 64, 0},
        };
        int gid = -1;
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < lanes.size(); ++i) {
            const auto& pl = lanes[i];
            std::vector<std::string> attrs = {
                pl.id, std::to_string((int)pl.turn_type), pl.entry_lane_id, pl.exit_lane_id,
                std::to_string((int)pl.lane_type),(pl.fixed_shape?"1":"0")
            };
            // points
            std::vector<ShapePoint> points;
            if (!pl.geometry.points.empty()) {
                points = toShapepoints(pl.geometry.points);
            } else {
                if (!pl.curve || pl.curve->numSegments() == 0)
                    continue;
                points = toShapepoints(pl.curve->sampleByShape());
            }
            ShapeRecord record = {
                ++gid, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeLaneEdges = [&](const std::string& dir, const std::string& fname,
                              const std::vector<ConnectivityLaneEdge>& edgelines) {
        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"GROUP_ID", 'C', 64, 0},
            {"GROUP_TYPE", 'C', 64, 0},
            {"LEFT_CLINE_ID", 'C', 64, 0},
            {"RIGHT_CLINE_ID", 'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < edgelines.size(); ++i) {
            const auto& pl = edgelines[i];
            std::vector<std::string> attrs = {pl.id, "-1", "", "", ""};
            std::vector<ShapePoint> points = toShapepoints(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeAreas = [&](const std::string& dir, const std::string& fname,
                          const IntersectionArea& inter_area) {
        std::vector<Polygon2d> areas = {inter_area.geometry};

        fs::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID", 'C', 64, 0},
            {"TYPE", 'C', 64, 0},
        };
        int shpType = SHP_POLYGONZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < areas.size(); ++i) {
            const auto& pg = areas[i];
            std::vector<std::string> attrs = {std::to_string(i), "-1"}; //std::to_string(pg.type)
            std::vector<ShapePoint> points = toShapepoints(pg.outer);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
#ifdef PROJECT_ROOT_DIR
    std::string proj_dir = PROJECT_ROOT_DIR;
    std::string qgis_fname = "intersection_gen.qgs";
    std::string qgis_temp = proj_dir + "/datas/" + qgis_fname;
    if (!fs::exists(qgis_temp)) qgis_temp = proj_dir + "/" + qgis_fname;
    if (fs::exists(qgis_temp)) {
        std::string dest_file = fs::path(out_dir).parent_path().parent_path().string() + "/" +qgis_fname;
        fs::copy_file(qgis_temp, dest_file, fs::copy_options::skip_existing);
    }
#endif
    writeLanes(out_dir, prefix + "_lanes", out.connectivity_curves);
    writeLaneEdges(out_dir, prefix + "_laneedges", out.lane_edges);
    writeAreas(out_dir, prefix + "_areas", out.area);
    return true;
}

}
#endif //IODATA_SHAPEFILE_H
