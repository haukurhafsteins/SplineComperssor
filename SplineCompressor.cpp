#pragma once
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <functional>
#include <vector>
#include <algorithm>
#include "SplineCompressor.hpp"

float SplineCompressor::estimateEpsilon(const float* array, size_t size) {
    if (size < 2) return 0.1f;

    float minVal = array[0], maxVal = array[0], sum = 0.0f, sumSq = 0.0f;

    for (size_t i = 0; i < size; i++) {
        float v = array[i];
        if (v < minVal) minVal = v;
        if (v > maxVal) maxVal = v;
        sum += v;
        sumSq += v * v;
    }

    float range = maxVal - minVal;
    float mean = sum / size;
    float variance = sumSq / size - mean * mean;
    float stdDev = std::sqrt(variance);

    float epsilon = std::fmax(0.05f * range, 0.5f * stdDev);
    if (epsilon < 0.1f) epsilon = 0.1f;
    if (epsilon > 0.25f * range) epsilon = 0.25f * range;

    return epsilon;
}

size_t SplineCompressor::rdpSimplify(const float* array, size_t size, Point* out, float epsilon, size_t max_out_size) {
    if (size < 2 || !out || max_out_size < 2) return 0;

    std::vector<bool> keep(size, false);
    std::vector<Point> points(size);

    for (size_t i = 0; i < size; i++) {
        points[i] = { static_cast<float>(i), array[i] };
    }

    keep[0] = true;
    keep[size - 1] = true;
    recursiveRdp(points, 0, size - 1, epsilon, keep);

    size_t count = 0;
    for (size_t i = 0; i < size && count < max_out_size; i++) {
        if (keep[i]) {
            out[count++] = points[i];
        }
    }

    return count;
}

size_t SplineCompressor::sampleSpline(const Point* control, size_t control_size,
    float* out, size_t max_out_points,
    float xStart, float xEnd) {
    if (control_size < 2 || max_out_points < 2 || !out) return 0;

    size_t n = control_size - 1;

    std::vector<float> a(control_size);
    std::vector<float> b(n);
    std::vector<float> c(control_size);
    std::vector<float> d(n);
    std::vector<float> h(n);
    std::vector<float> alpha(n);
    std::vector<float> l(control_size);
    std::vector<float> mu(control_size);
    std::vector<float> z(control_size);

    for (size_t i = 0; i < control_size; i++) a[i] = control[i].y;
    for (size_t i = 0; i < n; i++) h[i] = control[i + 1].x - control[i].x;

    for (size_t i = 1; i < n; i++) {
        alpha[i] = (3.0f / h[i]) * (a[i + 1] - a[i]) - (3.0f / h[i - 1]) * (a[i] - a[i - 1]);
    }

    l[0] = 1.0f; mu[0] = z[0] = 0.0f;

    for (size_t i = 1; i < n; i++) {
        l[i] = 2.0f * (h[i - 1] + h[i]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1.0f; z[n] = c[n] = 0.0f;

    for (int j = static_cast<int>(n) - 1; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0f * c[j]) / 3.0f;
        d[j] = (c[j + 1] - c[j]) / (3.0f * h[j]);
    }

    auto evalSpline = [&](size_t idx, float x) -> float {
        float dx = x - control[idx].x;
        return a[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx;
    };

    float minX = control[0].x;
    float maxX = control[n].x;
    float scale = (xEnd - xStart) / (maxX - minX);  // mapping factor

    float minY = std::numeric_limits<float>::max();
    float maxY = std::numeric_limits<float>::lowest();
    for (size_t i = 0; i < control_size; i++) {
        minY = std::min(minY, control[i].y);
        maxY = std::max(maxY, control[i].y);
    }
    float tolerance = std::max(0.01f * (maxY - minY), 0.001f); // adaptive density target

    auto mapX = [&](float x) -> float {
        return xStart + (x - minX) * scale;
    };

    std::vector<Point> sampled;
    sampled.reserve(max_out_points);
    sampled.push_back({ mapX(minX), control[0].y });

    std::function<void(size_t, float, float, float, float)> appendAdaptive =
        [&](size_t idx, float x0, float y0, float x1, float y1) {
            float midX = 0.5f * (x0 + x1);
            float midY = evalSpline(idx, midX);
            float linMid = y0 + (y1 - y0) * 0.5f;
            float error = std::fabs(midY - linMid);

            if (error > tolerance && sampled.size() + 1 < max_out_points) {
                appendAdaptive(idx, x0, y0, midX, midY);
                appendAdaptive(idx, midX, midY, x1, y1);
            } else {
                if (sampled.size() < max_out_points) {
                    sampled.push_back({ mapX(x1), y1 });
                } else {
                    sampled.back() = { mapX(x1), y1 };
                }
            }
        };

    for (size_t i = 0; i < n && sampled.size() < max_out_points; i++) {
        float x0 = control[i].x;
        float y0 = control[i].y;
        float x1 = control[i + 1].x;
        float y1 = evalSpline(i, x1); // ensure continuity with spline coefficients
        appendAdaptive(i, x0, y0, x1, y1);
    }

    size_t outCount = std::min(sampled.size(), max_out_points);
    for (size_t i = 0; i < outCount; i++) {
        out[i * 2] = sampled[i].x;
        out[i * 2 + 1] = sampled[i].y;
    }

    return outCount;
}

//---------------------------------------------------------
// Private functions
//---------------------------------------------------------
void SplineCompressor::recursiveRdp(std::vector<Point>& pts, size_t start, size_t end, float epsilon, std::vector<bool>& keep) {
    if (end <= start + 1) return;

    float maxDist = 0.0f;
    size_t index = start;
    for (size_t i = start + 1; i < end; i++) {
        float dist = perpendicularDistance(pts[i], pts[start], pts[end]);
        if (dist > maxDist) {
            maxDist = dist;
            index = i;
        }
    }

    if (maxDist > epsilon) {
        keep[index] = true;
        recursiveRdp(pts, start, index, epsilon, keep);
        recursiveRdp(pts, index, end, epsilon, keep);
    }
}

float SplineCompressor::perpendicularDistance(Point p, Point start, Point end) {
    float dx = end.x - start.x;
    float dy = end.y - start.y;
    float num = std::fabs(dy * p.x - dx * p.y + end.x * start.y - end.y * start.x);
    float denom = std::sqrt(dx * dx + dy * dy);
    return num / (denom + 1e-6f);
}
