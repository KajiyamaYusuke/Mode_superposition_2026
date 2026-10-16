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
    if (count < 3 || valid.size() != areaMm2.size()) return result;

    const int lastInteriorIndex = count - 2;
    double minimum = std::numeric_limits<double>::infinity();
    for (int i = 1; i <= lastInteriorIndex; ++i) {
        if (!valid[i] || !std::isfinite(areaMm2[i])) continue;
        if (areaMm2[i] < minimum) {
            minimum = areaMm2[i];
            result.minIndex = i;
        }
    }
    if (result.minIndex < 0) return result;

    result.minAreaMm2 = minimum;
    result.targetAreaMm2 = areaRatio * minimum;

    if (minimum <= closedAreaMm2) {
        // Preserve the former closure convention: the first non-positive
        // valid section takes precedence; otherwise use the global minimum.
        result.sepIndex = result.minIndex;
        for (int i = 1; i <= lastInteriorIndex; ++i) {
            if (valid[i] && std::isfinite(areaMm2[i]) && areaMm2[i] <= 0.0) {
                result.sepIndex = i;
                break;
            }
        }
        result.sepAreaMm2 = areaMm2[result.sepIndex];
        result.reason = "closed_section";
        return result;
    }

    int lastValidInterior = result.minIndex;
    for (int i = result.minIndex; i <= lastInteriorIndex; ++i) {
        if (!valid[i] || !std::isfinite(areaMm2[i])) continue;
        lastValidInterior = i;
        if (areaMm2[i] >= result.targetAreaMm2) {
            result.sepIndex = i;
            result.sepAreaMm2 = areaMm2[i];
            result.reason = "area_ratio_reached";
            return result;
        }
    }

    result.sepIndex = lastValidInterior;
    result.sepAreaMm2 = areaMm2[lastValidInterior];
    result.usedFallback = true;
    result.reason = "no_downstream_ratio_crossing";
    return result;
}
