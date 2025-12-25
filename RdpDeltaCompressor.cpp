#include "RdpDeltaCompressor.hpp"
#include <cmath>
#include <limits>

size_t RdpDeltaCompressor::compress(const Point* points, size_t size, float epsilon, float scale, Point& first, DeltaPair* out, size_t max_out) {
    if (!points || size < 2 || !out || max_out == 0 || scale <= 0.0f) return 0;

    std::vector<Point> simplified = rdpSimplify(points, size, epsilon);
    if (simplified.size() < 2) return 0;
    if (simplified.size() - 1 > max_out) return 0;

    first = simplified[0];
    Point prev = first;

    int idx = 0;
    for (size_t i = 1; i < simplified.size(); i++) {
        float dx = (simplified[i].x - prev.x) * scale;
        float dy = (simplified[i].y - prev.y) * scale;
        long lx = std::lround(dx);
        long ly = std::lround(dy);
        if (lx < std::numeric_limits<int16_t>::min() || lx > std::numeric_limits<int16_t>::max()) return 0;
        if (ly < std::numeric_limits<int16_t>::min() || ly > std::numeric_limits<int16_t>::max()) return 0;
        out[idx].dx = static_cast<int16_t>(lx);
        out[idx].dy = static_cast<int16_t>(ly);
        idx++;
        prev = simplified[i];
    }

    return simplified.size();
}

size_t RdpDeltaCompressor::decompress(const Point& first, const DeltaPair* deltas, size_t delta_count, float scale, Point* out, size_t max_out) {
    if (!deltas || !out || delta_count + 1 > max_out || scale <= 0.0f) return 0;
    out[0] = first;
    for (size_t i = 0; i < delta_count; i++) {
        Point p;
        p.x = out[i].x + static_cast<float>(deltas[i].dx) / scale;
        p.y = out[i].y + static_cast<float>(deltas[i].dy) / scale;
        out[i + 1] = p;
    }
    return delta_count + 1;
}

std::vector<Point> RdpDeltaCompressor::rdpSimplify(const Point* points, size_t size, float epsilon) {
    std::vector<bool> keep(size, false);
    keep[0] = true;
    keep[size - 1] = true;

    std::vector<std::pair<size_t, size_t>> stack;
    stack.reserve(size);
    stack.push_back({0, size - 1});

    while (!stack.empty()) {
        auto [start, end] = stack.back();
        stack.pop_back();

        float maxDist = 0.0f;
        size_t index = start;
        for (size_t i = start + 1; i < end; i++) {
            float dist = perpendicularDistance(points[i], points[start], points[end]);
            if (dist > maxDist) {
                maxDist = dist;
                index = i;
            }
        }

        if (maxDist > epsilon) {
            keep[index] = true;
            stack.push_back({start, index});
            stack.push_back({index, end});
        }
    }

    std::vector<Point> simplified;
    simplified.reserve(size);
    for (size_t i = 0; i < size; i++) {
        if (keep[i]) simplified.push_back(points[i]);
    }
    return simplified;
}

float RdpDeltaCompressor::perpendicularDistance(Point p, Point start, Point end) {
    float dx = end.x - start.x;
    float dy = end.y - start.y;
    float num = std::fabs(dy * p.x - dx * p.y + end.x * start.y - end.y * start.x);
    float denom = std::sqrt(dx * dx + dy * dy);
    return num / (denom + 1e-6f);
}
