#include "ChannelSections.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

constexpr double intersectionTolerance = 1.0e-12;
struct SectionPoint { double x, y, z; };
struct SectionSegment { SectionPoint a, b; };
struct GapIntegrals { double area = 0.0, cubed = 0.0, span = 0.0; };

std::pair<double, double> surfaceXRange(const Geometry& geom, const State& state) {
    if (geom.nxsup < 2 || geom.nxsup > geom.nsurfl || geom.nsurfz < 2
        || static_cast<int>(geom.surfp.size()) < geom.nxsup)
        throw std::runtime_error("invalid channel surface dimensions");
    double low = std::numeric_limits<double>::infinity(), high = -low;
    for (int i = 0; i < geom.nxsup; ++i) {
        if (static_cast<int>(geom.surfp[i].size()) < geom.nsurfz)
            throw std::runtime_error("incomplete channel surface row");
        for (int j = 0; j < geom.nsurfz; ++j) {
            const int id = geom.surfp[i][j];
            if (id == -1) continue;
            if (id < 0 || id >= static_cast<int>(state.disp.size()))
                throw std::runtime_error("invalid channel surface node index");
            const auto& p = state.disp[id]; // Absolute deformed coordinates.
            if (!std::isfinite(p.ux) || !std::isfinite(p.uy) || !std::isfinite(p.uz))
                throw std::runtime_error("nonfinite channel surface node " + std::to_string(id));
            low = std::min(low, p.ux); high = std::max(high, p.ux);
        }
    }
    return {low, high};
}

bool intersectAtX(const Geometry& geom, const State& state, int j,
                  double xPlane, bool takeMaximumY, SectionPoint& result) {
    std::vector<SectionPoint> candidates;
    auto consider = [&](SectionPoint p) {
        for (const auto& q : candidates) {
            if (std::hypot(p.y - q.y, p.z - q.z) <= intersectionTolerance) return;
        }
        candidates.push_back(p);
        if (candidates.size() == 1 || (takeMaximumY ? p.y > result.y : p.y < result.y))
            result = p;
    };
    for (int i = 0; i + 1 < geom.nxsup; ++i) {
        const int ia = geom.surfp[i][j], ib = geom.surfp[i + 1][j];
        if (ia < 0 || ib < 0) continue; // Never bridge a missing surface node.
        const auto& a = state.disp[ia]; const auto& b = state.disp[ib];
        if (xPlane < std::min(a.ux, b.ux) - intersectionTolerance
            || xPlane > std::max(a.ux, b.ux) + intersectionTolerance) continue;
        const double dx = b.ux - a.ux;
        if (std::abs(dx) <= intersectionTolerance) {
            // A plane-aligned edge contributes its inner endpoint, not an
            // arbitrary midpoint or the segment nearest a predicted i.
            if (std::abs(xPlane - a.ux) <= intersectionTolerance
                && std::abs(xPlane - b.ux) <= intersectionTolerance) {
                consider({xPlane, a.uy, a.uz});
                consider({xPlane, b.uy, b.uz});
            }
        } else {
            const double s = std::clamp((xPlane - a.ux) / dx, 0.0, 1.0);
            consider({xPlane, a.uy + s * (b.uy - a.uy), a.uz + s * (b.uz - a.uz)});
        }
    }
    return !candidates.empty();
}

std::vector<SectionSegment> sectionSegments(const Geometry& geom, const State& state,
                                           double x, bool takeMaximumY, const char* side) {
    std::vector<SectionSegment> segments;
    SectionPoint previous{}, lastFound{};
    bool previousValid = false, haveFound = false;
    for (int j = 0; j < geom.nsurfz; ++j) {
        SectionPoint p{};
        const bool found = intersectAtX(geom, state, j, x, takeMaximumY, p);
        if (found) {
            if (haveFound && p.z < lastFound.z - intersectionTolerance)
                throw std::runtime_error(std::string(side) + " AP z order reversal at j="
                    + std::to_string(j) + " x_mm=" + std::to_string(x));
            // Only consecutive valid j values define an interpolatable AP
            // segment. Horizontal boundaries have zero integration width.
            if (previousValid && p.z > previous.z) segments.push_back({previous, p});
            previous = p; lastFound = p; haveFound = true;
        }
        previousValid = found;
    }
    return segments;
}

double interpolateYAtZ(const SectionSegment& segment, double z) {
    const double s = (z - segment.a.z) / (segment.b.z - segment.a.z);
    return segment.a.y + s * (segment.b.y - segment.a.y);
}

void positiveLinearIntegrals(double a, double b, double length, double& area, double& cubed) {
    // Keep signed endpoints until the exact zero crossing has been split:
    // negative gap is contact, never an invalid geometric search.
    if (a <= 0.0 && b <= 0.0) return;
    if (a < 0.0) { length *= b / (b - a); a = 0.0; }
    else if (b < 0.0) { length *= a / (a - b); b = 0.0; }
    area += 0.5 * (a + b) * length;
    cubed += 0.25 * (a*a*a + a*a*b + a*b*b + b*b*b) * length;
}


GapIntegrals positiveGapIntegralsYZ(const std::vector<SectionSegment>& left,
                                    const std::vector<SectionSegment>& right,
                                    double gapSign) {
    GapIntegrals out;
    size_t il = 0, ir = 0;
    // Ordered segment intersections split at every actual left/right z
    // vertex, and sum disjoint common AP intervals without bridging holes.
    while (il < left.size() && ir < right.size()) {
        const auto& l = left[il]; const auto& r = right[ir];
        const double low = std::max(l.a.z, r.a.z), high = std::min(l.b.z, r.b.z);
        if (high > low) {
            const double h0 = gapSign * (interpolateYAtZ(r, low) - interpolateYAtZ(l, low));
            const double h1 = gapSign * (interpolateYAtZ(r, high) - interpolateYAtZ(l, high));
            if (!std::isfinite(h0) || !std::isfinite(h1))
                throw std::runtime_error("nonfinite signed channel gap");
            positiveLinearIntegrals(h0, h1, high - low, out.area, out.cubed);
            out.span += high - low;
        }
        if (l.b.z <= r.b.z) ++il;
        if (r.b.z <= l.b.z) ++ir;
    }
    return out;
}
} // namespace

void ChannelSections::resize(const Geometry& left, const Geometry& right, int nSections) {
    if (nSections < 3) throw std::runtime_error("channel requires at least three sections");
    // Compatibility fields only; no unused material anchor affects the x grid.
    sigma.clear();
    referenceConstrictionI = -1;
    referenceConstrictionSigma = std::numeric_limits<double>::quiet_NaN();
    x.assign(nSections, 0.0);
    area.assign(nSections, std::numeric_limits<double>::quiet_NaN());
    gapCubed = area;
    valid.assign(nSections, false);
    double meanInitialGap = 0.0;
    int count = 0;
    for (int i = 0; i < std::min(left.nxsup, right.nxsup); ++i) {
        for (int j = 0; j < std::min(left.nsurfz, right.nsurfz); ++j) {
            const int il = left.surfp.at(i).at(j), ir = right.surfp.at(i).at(j);
            if (il < 0 || ir < 0) continue;
            meanInitialGap += right.points.at(ir).y - left.points.at(il).y;
            ++count;
        }
    }
    if (count == 0 || !std::isfinite(meanInitialGap))
        throw std::runtime_error("cannot determine initial channel gap orientation");
    gapSign = meanInitialGap < 0.0 ? -1.0 : 1.0;
}

void ChannelSectionBuilder::build(const Geometry& leftGeom, const State& leftState,
                                  const Geometry& rightGeom, const State& rightState,
                                  ChannelSections& sections) const {
    const size_t count = sections.x.size();
    if (count < 3 || sections.area.size() != count || sections.gapCubed.size() != count
        || sections.valid.size() != count) throw std::runtime_error("inconsistent channel arrays");
    const auto l = surfaceXRange(leftGeom, leftState);
    const auto r = surfaceXRange(rightGeom, rightState);
    // Intersect the two WHOLE surface x ranges, not all individual j ranges.
    const double low = std::max(l.first, r.first), high = std::min(l.second, r.second);
    if (!std::isfinite(low) || !std::isfinite(high) || !(high > low))
        throw std::runtime_error("no finite common channel x range");
    std::fill(sections.valid.begin(), sections.valid.end(), false);
    std::fill(sections.area.begin(), sections.area.end(), std::numeric_limits<double>::quiet_NaN());
    std::fill(sections.gapCubed.begin(), sections.gapCubed.end(), std::numeric_limits<double>::quiet_NaN());
    for (size_t k = 0; k < count; ++k) {
        sections.x[k] = low + (high - low) * static_cast<double>(k) / (count - 1);
        if (!std::isfinite(sections.x[k]) || (k > 0 && !(sections.x[k] > sections.x[k - 1])))
            throw std::runtime_error("nonmonotonic channel x section " + std::to_string(k));
        const auto left = sectionSegments(leftGeom, leftState, sections.x[k], sections.gapSign > 0.0, "left");
        const auto right = sectionSegments(rightGeom, rightState, sections.x[k], sections.gapSign < 0.0, "right");
        const auto integral = positiveGapIntegralsYZ(left, right, sections.gapSign);
        if (!std::isfinite(integral.area) || !std::isfinite(integral.cubed)
            || integral.area < 0.0 || integral.cubed < 0.0)
            throw std::runtime_error("invalid channel integral at section " + std::to_string(k));
        if (integral.span > 0.0) {
            sections.area[k] = integral.area;
            sections.gapCubed[k] = integral.cubed;
            sections.valid[k] = true;
        }
    }
}
