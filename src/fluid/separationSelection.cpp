#include "SeparationSelection.h"

#include <cmath>
#include <limits>

SeparationSelection selectSeparationByAreaRatio(
    const std::vector<double>& areaMm2,
    const std::vector<bool>& valid,
    double areaRatio,
    double closedAreaMm2) {
    SeparationSelection result;
    const int count = static_cast<int>(areaMm2.size());
    if (count < 2 || valid.size() != areaMm2.size()) return result;
    if (!std::isfinite(areaRatio) || areaRatio < 1.0
        || !std::isfinite(closedAreaMm2) || closedAreaMm2 < 0.0) {
        result.reason = "invalid_separation_parameters";
        return result;
    }

    // Missing geometry ends the domain. Never restart at a later valid plane.
    int end = 0;
    while (end < count && valid[end]) {
        if (!std::isfinite(areaMm2[end]) || areaMm2[end] < 0.0) {
            result.reason = "invalid_area_section_" + std::to_string(end);
            return result;
        }
        ++end;
    }
    if (end < 2) {
        result.reason = "insufficient_inlet_coverage";
        return result;
    }
    const int last = end - 1;
    const bool truncated = end < count;
    double minimum = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= last; ++i) {
        if (areaMm2[i] < minimum) {
            minimum = areaMm2[i];
            result.minIndex = i;
        }
    }
    result.minAreaMm2 = minimum;
    result.targetAreaMm2 = areaRatio * minimum;
    if (!std::isfinite(result.targetAreaMm2)) {
        result.reason = "nonfinite_target_area";
        return result;
    }
    if (truncated && result.minIndex == last) {
        result.reason = "coverage_limited_minimum_section_" + std::to_string(last);
        return result;
    }

    if (minimum <= closedAreaMm2) {
        result.sepIndex = result.minIndex;
        // Retain the true-closure convention, without treating missing cuts
        // as zero area. Negative signed gaps were clipped only in A and G3.
        for (int i = 0; i <= last; ++i) {
            if (areaMm2[i] <= 0.0) { result.sepIndex = i; break; }
        }
        result.sepAreaMm2 = areaMm2[result.sepIndex];
        result.reason = "closed_section";
        return result;
    }

    for (int i = result.minIndex; i <= last; ++i) {
        if (areaMm2[i] >= result.targetAreaMm2) {
            if (truncated && i == last) {
                result.reason = "coverage_limited_separation_section_" + std::to_string(last);
                return result;
            }
            result.sepIndex = i;
            result.sepAreaMm2 = areaMm2[i];
            result.reason = "area_ratio_reached";
            return result;
        }
    }
    if (truncated) {
        result.reason = "coverage_limited_no_ratio_crossing_section_" + std::to_string(last);
        return result;
    }
    // Only the complete configured domain permits a last-plane fallback.
    result.sepIndex = last;
    result.sepAreaMm2 = areaMm2[last];
    result.usedFallback = true;
    result.reason = "no_downstream_ratio_crossing";
    return result;
}
