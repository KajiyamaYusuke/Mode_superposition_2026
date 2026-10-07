#include "SeparationSelection.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
int failures = 0;

void expect(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << "\n";
        ++failures;
    }
}

void expectNear(double actual, double expected, const std::string& message) {
    expect(std::abs(actual - expected) < 1.0e-12, message);
}
}

int main() {
    constexpr double closed = 0.01;

    {
        const std::vector<double> area{4.0, 2.0, 1.0, 1.1, 1.19, 1.21, 1.5};
        const auto result = selectSeparationByAreaRatio(
            area, std::vector<bool>(area.size(), true), 1.2, closed);
        expect(result.minIndex == 2, "normal: minIndex");
        expect(result.sepIndex == 5, "normal: sepIndex");
        expectNear(result.targetAreaMm2, 1.2, "normal: target area");
        expect(!result.usedFallback, "normal: no fallback");
        expect(result.reason == "area_ratio_reached", "normal: reason");
    }
    {
        const std::vector<double> area{4.0, 2.0, 1.0, 1.2, 1.4};
        const auto result = selectSeparationByAreaRatio(
            area, std::vector<bool>(area.size(), true), 1.2, closed);
        expect(result.sepIndex == 3, "exact threshold");
    }
    {
        const std::vector<double> area{4.0, 2.0, 1.0, 1.1, 1.4};
        const auto result = selectSeparationByAreaRatio(
            area, std::vector<bool>(area.size(), true), 1.0, closed);
        expect(result.sepIndex == result.minIndex, "ratio 1.0 selects minimum");
    }
    {
        const std::vector<double> area{4.0, 2.0, 1.0, 1.05, 1.1};
        const auto result = selectSeparationByAreaRatio(
            area, std::vector<bool>(area.size(), true), 1.2, closed);
        expect(result.sepIndex == 3, "fallback selects last valid interior");
        expect(result.usedFallback, "fallback flag");
        expect(result.reason == "no_downstream_ratio_crossing", "fallback reason");
    }
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const std::vector<double> area{4.0, 2.0, 1.0, 1.25, nan, 1.3, 1.5};
        std::vector<bool> valid(area.size(), true);
        valid[3] = false;
        const auto result = selectSeparationByAreaRatio(area, valid, 1.2, closed);
        expect(result.sepIndex == 5, "invalid and NaN sections are skipped");
    }
    {
        const std::vector<double> area{4.0, 0.005, 0.0, 0.004, 1.0};
        const auto result = selectSeparationByAreaRatio(
            area, std::vector<bool>(area.size(), true), 1.2, closed);
        expect(result.minIndex == 2, "closure: global minimum");
        expect(result.sepIndex == 2, "closure: no ratio search");
        expect(!result.usedFallback, "closure: no fallback");
        expect(result.reason == "closed_section", "closure reason");
    }

    return failures == 0 ? 0 : 1;
}
