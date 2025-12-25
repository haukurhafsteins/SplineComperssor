#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include "SplineCompressor.hpp"

struct DeltaPair {
    int16_t dx;
    int16_t dy;
};

class RdpDeltaCompressor {
public:
    /// @brief Simplify a point series with RDP and delta-encode the result.
    /// @param points   Input points.
    /// @param size     Number of input points.
    /// @param epsilon  RDP tolerance.
    /// @param scale    Fixed-point scale (value * scale -> stored as int16).
    /// @param first    Output: anchor point stored in full precision.
    /// @param out      Output: delta pairs (dx, dy) after the anchor.
    /// @param max_out  Capacity of @p out (number of delta pairs, not floats).
    /// @return Number of kept points (anchor + deltas), or 0 on failure.
    static size_t compress(const Point* points, size_t size, float epsilon, float scale, Point& first, DeltaPair* out, size_t max_out);

    /// @brief Reconstruct the simplified points from the anchor and delta pairs.
    /// @param first        Anchor point.
    /// @param deltas       Delta pairs.
    /// @param delta_count  Number of delta pairs.
    /// @param scale        Fixed-point scale used during compression.
    /// @param out          Output buffer for reconstructed points.
    /// @param max_out      Capacity of @p out (points).
    /// @return Number of reconstructed points, or 0 on failure.
    static size_t decompress(const Point& first, const DeltaPair* deltas, size_t delta_count, float scale, Point* out, size_t max_out);

private:
    static std::vector<Point> rdpSimplify(const Point* points, size_t size, float epsilon);
    static float perpendicularDistance(Point p, Point start, Point end);
};
