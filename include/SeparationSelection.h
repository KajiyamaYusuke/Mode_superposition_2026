#pragma once

#include <vector>
#include <string>

struct SeparationSelection {
    int minIndex = -1;
    int sepIndex = -1;
    double minAreaMm2 = 0.0;
    double targetAreaMm2 = 0.0;
    double sepAreaMm2 = 0.0;
    bool usedFallback = false;
    std::string reason = "invalid_sections";
};

SeparationSelection selectSeparationByAreaRatio(
    const std::vector<double>& areaMm2,
    const std::vector<bool>& valid,
    double areaRatio,
    double closedAreaMm2);
