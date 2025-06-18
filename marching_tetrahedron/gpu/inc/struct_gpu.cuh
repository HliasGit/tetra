#pragma once

#include <cuda_runtime.h>
#include <unordered_set>
#include "global.h"

typedef struct cube_gpu_SoA{
    float4 *coord_idx;
    float4 *one_apex;
    float4 *two_apex;
} cube_gpu_SoA;

typedef struct cube_vertices_points_SoA{
    float4 val;
} cube_vertices_points_SoA;

typedef struct TriangleVertex_GPU{
    coord_t x;
    coord_t y;
    coord_t z;

    bool operator==(const TriangleVertex_GPU& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator<(const TriangleVertex_GPU& other) const{
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
} TriangleVertex_GPU;

typedef struct Triangle_GPU{
    TriangleVertex_GPU v1;
    TriangleVertex_GPU v2;
    TriangleVertex_GPU v3;
} Triangle_GPU;

namespace std {
    template <>
    struct hash<TriangleVertex_GPU> {
        std::size_t operator()(const TriangleVertex_GPU& v) const {
            std::size_t hx = std::hash<coord_t>{}(v.x);
            std::size_t hy = std::hash<coord_t>{}(v.y);
            std::size_t hz = std::hash<coord_t>{}(v.z);
            return hx ^ (hy << 1) ^ (hz << 2);
        }
    };
}
