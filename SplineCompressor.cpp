#pragma once
#include <cstdio>
#include <cmath>
#include <cstddef>
#include <cstdlib>
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

    bool* keep = (bool*)calloc(size, sizeof(bool));
    Point* points = (Point*)malloc(size * sizeof(Point));
    if (!keep || !points) {
        free(keep); free(points);
        return 0;
    }

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

    free(keep);
    free(points);
    return count;
}

bool SplineCompressor::sampleSpline(const Point* control, size_t control_size,
    float* out, size_t sample_count,
    float xStart, float xEnd) {
    if (control_size < 2 || sample_count < 2 || !out) return false;

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

    float minX = control[0].x;
    float maxX = control[n].x;
    float step = (maxX - minX) / (sample_count - 1);
    float scale = (xEnd - xStart) / (maxX - minX);  // mapping factor

    for (size_t i = 0; i < sample_count; i++) {
        float x = minX + i * step;
        size_t idx = findSegment(control, n, x);
        float dx = x - control[idx].x;

        float y = a[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx;

        // Map x to output domain
        out[i * 2] = xStart + (x - minX) * scale;
        out[i * 2 + 1] = y;
    }

    return true;
}

//---------------------------------------------------------
// Private functions
//---------------------------------------------------------
void SplineCompressor::recursiveRdp(Point* pts, size_t start, size_t end, float epsilon, bool* keep) {
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

size_t SplineCompressor::findSegment(const Point* points, size_t n, float x) {
    for (size_t i = 0; i < n; i++) {
        if (x < points[i + 1].x) return i;
    }
    return n - 1;
}
