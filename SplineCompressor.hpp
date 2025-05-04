#pragma once
#include <cmath>
#include <cstddef>
#include <cstdlib>

struct Point {
    float x;
    float y;
};

class SplineCompressor {
public:
    static float estimateEpsilon(const float* array, size_t size);
    static size_t rdpSimplify(const float* array, size_t size, Point* out, float epsilon, size_t max_out_size);
    static bool sampleSpline(const Point* control, size_t control_size, float* out, size_t sample_count, float xStart, float xEnd);

private:
    static void recursiveRdp(Point* pts, size_t start, size_t end, float epsilon, bool* keep);
    static float perpendicularDistance(Point p, Point start, Point end);

    static size_t findSegment(const Point* points, size_t n, float x);
};
