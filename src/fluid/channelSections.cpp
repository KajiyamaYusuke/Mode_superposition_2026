#include "ChannelSections.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

struct SectionPoint { double x, y, z; };

bool intersectAtX(const Geometry& geom, const State& state, int j,
                  double xPlane, int preferredSegment, SectionPoint& result) {
    const int n = geom.nxsup;
    if (n < 2 || j < 0 || j >= geom.nsurfz) return false;

    int best = -1;
    double bestDistance = std::numeric_limits<double>::infinity();
    for (int k = 0; k < n - 1; ++k) {
        const int a = geom.surfp[k][j];
        const int b = geom.surfp[k + 1][j];
        if (a < 0 || b < 0) continue;
        const double xa = state.disp[a].ux;
        const double xb = state.disp[b].ux;
        if ((xPlane - xa) * (xPlane - xb) > 0.0) continue;
        const double distance = std::abs(k - preferredSegment);
        if (distance < bestDistance) { best = k; bestDistance = distance; }
    }
    if (best < 0) return false;

    const int a = geom.surfp[best][j];
    const int b = geom.surfp[best + 1][j];
    const auto& pa = state.disp[a];
    const auto& pb = state.disp[b];
    const double dx = pb.ux - pa.ux;
    const double t = std::abs(dx) > 1e-12 ? (xPlane - pa.ux) / dx : 0.5;
    result = {xPlane, pa.uy + t * (pb.uy - pa.uy), pa.uz + t * (pb.uz - pa.uz)};
    return std::isfinite(result.y) && std::isfinite(result.z);
}

struct GapIntegrals { double area = 0.0, cubed = 0.0; };

double interpolateYAtZ(const std::vector<SectionPoint>& curve, double z) {
    for (size_t i = 0; i + 1 < curve.size(); ++i) {
        const double z0 = curve[i].z, z1 = curve[i + 1].z;
        if ((z - z0) * (z - z1) > 0.0) continue;
        const double dz = z1 - z0;
        const double s = std::abs(dz) > 1.0e-12 ? (z - z0) / dz : 0.5;
        return curve[i].y + s * (curve[i + 1].y - curve[i].y);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

GapIntegrals positiveGapIntegralsYZ(std::vector<SectionPoint> left,
                                    std::vector<SectionPoint> right,
                                    double gapSign) {
    if (left.size() < 2 || right.size() < 2) return {};
    auto byZ = [](const SectionPoint& a, const SectionPoint& b) { return a.z < b.z; };
    std::sort(left.begin(), left.end(), byZ);
    std::sort(right.begin(), right.end(), byZ);
    const double zLow = std::max(left.front().z, right.front().z);
    const double zHigh = std::min(left.back().z, right.back().z);
    if (!(zHigh > zLow)) return {};

    // A shared physical z grid prevents a deformed left j-line from being
    // integrated against a different physical point on the right surface.
    std::vector<double> zGrid{zLow, zHigh};
    for (const auto& p : left)  if (p.z > zLow && p.z < zHigh) zGrid.push_back(p.z);
    for (const auto& p : right) if (p.z > zLow && p.z < zHigh) zGrid.push_back(p.z);
    std::sort(zGrid.begin(), zGrid.end());
    zGrid.erase(std::unique(zGrid.begin(), zGrid.end(), [](double a, double b) {
        return std::abs(a - b) < 1.0e-12;
    }), zGrid.end());

    GapIntegrals out;
    for (size_t k = 0; k + 1 < zGrid.size(); ++k) {
        const double z0 = zGrid[k], z1 = zGrid[k + 1], dz = z1 - z0;
        const double g0 = gapSign * (interpolateYAtZ(right, z0) - interpolateYAtZ(left, z0));
        const double g1 = gapSign * (interpolateYAtZ(right, z1) - interpolateYAtZ(left, z1));
        if (!(dz > 0.0) || !std::isfinite(g0) || !std::isfinite(g1)) {
            return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        }
        auto integratePositiveLinear = [&](double a, double b, double length) {
            GapIntegrals part;
            if (a >= 0.0 && b >= 0.0) {
                part.area = 0.5 * (a + b) * length;
                part.cubed = 0.25 * (a*a*a + a*a*b + a*b*b + b*b*b) * length;
            } else if (a >= 0.0 && b < 0.0) {
                const double positiveLength = length * a / (a - b);
                part.area = 0.5 * a * positiveLength;
                part.cubed = 0.25 * a*a*a * positiveLength;
            } else if (a < 0.0 && b >= 0.0) {
                const double positiveLength = length * b / (b - a);
                part.area = 0.5 * b * positiveLength;
                part.cubed = 0.25 * b*b*b * positiveLength;
            }
            return part;
        };
        const auto part = integratePositiveLinear(g0, g1, dz);
        out.area += part.area;
        out.cubed += part.cubed;
    }
    return out;
}
} // namespace

void ChannelSections::resize(const Geometry& left, const Geometry& right, int nSections) {
    const int nx = nSections;
    const int nzL = left.nsurfz;
    const int nzR = right.nsurfz;
    sigma.assign(nx, 0.0);
    x.assign(nx, 0.0);
    area.assign(nx, 0.0);
    gapCubed.assign(nx, 0.0);
    valid.assign(nx, false);
    for (int k = 0; k < nx; ++k) sigma[k] = nx > 1 ? static_cast<double>(k) / (nx - 1) : 0.0;

    // Preserve the narrowest initial i-grid slice as an explicit flow plane.
    // This is only used to choose the fixed material labels; all subsequent
    // areas are still calculated from interpolated intersections with the
    // deformed x = const planes below.
    const int nCommonI = std::min(left.nxsup, right.nxsup);
    const int nCommonJ = std::min(left.nsurfz, right.nsurfz);
    double meanInitialGap = 0.0;
    int initialGapCount = 0;
    for (int i = 0; i < nCommonI; ++i) {
        for (int j = 0; j < nCommonJ; ++j) {
            const int il = left.surfp[i][j];
            const int ir = right.surfp[i][j];
            if (il < 0 || ir < 0) continue;
            meanInitialGap += right.points[ir].y - left.points[il].y;
            ++initialGapCount;
        }
    }
    gapSign = (initialGapCount > 0 && meanInitialGap < 0.0) ? -1.0 : 1.0;

    double minimumArea = std::numeric_limits<double>::infinity();
    for (int i = 0; i < nCommonI; ++i) {
        std::vector<SectionPoint> leftSlice, rightSlice;
        leftSlice.reserve(nCommonJ);
        rightSlice.reserve(nCommonJ);
        bool complete = true;
        for (int j = 0; j < nCommonJ; ++j) {
            const int il = left.surfp[i][j];
            const int ir = right.surfp[i][j];
            if (il < 0 || ir < 0) { complete = false; break; }
            leftSlice.push_back({left.points[il].x, left.points[il].y, left.points[il].z});
            rightSlice.push_back({right.points[ir].x, right.points[ir].y, right.points[ir].z});
        }
        if (!complete) continue;
        const auto integrals = positiveGapIntegralsYZ(leftSlice, rightSlice, gapSign);
        if (std::isfinite(integrals.area) && integrals.area < minimumArea) {
            minimumArea = integrals.area;
            referenceConstrictionI = i;
        }
    }

    if (referenceConstrictionI >= 0 && nCommonI > 1 && nx > 2) {
        referenceConstrictionSigma =
            static_cast<double>(referenceConstrictionI) / (nCommonI - 1);

        // Replace the closest regular station, then sort so that all later
        // pressure interpolation remains ordered in material coordinate.
        const int kAnchor = std::clamp(
            static_cast<int>(std::lround(referenceConstrictionSigma * (nx - 1))),
            1, nx - 2);
        sigma[kAnchor] = referenceConstrictionSigma;
        std::sort(sigma.begin(), sigma.end());
    }
}

void ChannelSectionBuilder::build(const Geometry& leftGeom, const State& leftState,
                                  const Geometry& rightGeom, const State& rightState,
                                  ChannelSections& sections) const {
    const int nx = static_cast<int>(sections.sigma.size());
    const int nz = std::min(leftGeom.nsurfz, rightGeom.nsurfz);
    if (nx < 2 || nz < 2) return;

    // Every requested section is clamped to the x range common to all surface
    // flow lines.  A failed cut is therefore never silently converted to A=0.
    double commonLow = -std::numeric_limits<double>::infinity();
    double commonHigh = std::numeric_limits<double>::infinity();
    for (const auto& pair : {std::make_pair(&leftGeom, &leftState), std::make_pair(&rightGeom, &rightState)}) {
        const Geometry& geom = *pair.first;
        const State& state = *pair.second;
        for (int j = 0; j < nz; ++j) {
            double low = std::numeric_limits<double>::infinity();
            double high = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < geom.nxsup; ++i) {
                const int id = geom.surfp[i][j];
                if (id < 0) continue;
                low = std::min(low, state.disp[id].ux);
                high = std::max(high, state.disp[id].ux);
            }
            commonLow = std::max(commonLow, low);
            commonHigh = std::min(commonHigh, high);
        }
    }

    for (int i = 0; i < nx; ++i) {
        if (!(commonLow <= commonHigh) || !std::isfinite(commonLow) || !std::isfinite(commonHigh)) {
            sections.x[i] = std::numeric_limits<double>::quiet_NaN();
            sections.valid[i] = false;
            sections.area[i] = std::numeric_limits<double>::quiet_NaN();
            sections.gapCubed[i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        // Equal spacing in the instantaneous common physical range guarantees
        // strictly ordered, non-duplicated flow planes even when material
        // planes near an end clamp to the same x coordinate.
        sections.x[i] = commonLow + (commonHigh - commonLow)
            * (nx > 1 ? static_cast<double>(i) / (nx - 1) : 0.0);
        std::vector<SectionPoint> left, right;
        left.reserve(nz); right.reserve(nz);
        bool complete = true;
        const int nCommonI = std::min(leftGeom.nxsup, rightGeom.nxsup);
        const int preferredSegment = std::clamp(
            static_cast<int>(std::lround((sections.x[i] - commonLow)
                / std::max(commonHigh - commonLow, 1.0e-12) * (nCommonI - 1))),
            0, std::max(0, nCommonI - 2));
        for (int j = 0; j < nz; ++j) {
            SectionPoint pl{}, pr{};
            // The preferred segment must be derived from the material sigma,
            // not from the 50-section index.  Otherwise a 50-section flow
            // grid on a differently sized surface mesh changes branches when
            // a tilted surface intersects a plane more than once.
            complete = intersectAtX(leftGeom, leftState, j, sections.x[i], preferredSegment, pl)
                    && intersectAtX(rightGeom, rightState, j, sections.x[i], preferredSegment, pr);
            if (!complete) break;
            left.push_back(pl); right.push_back(pr);
        }
        sections.valid[i] = complete;
        // Invalid is distinct from a genuine closed section.  Keep a benign
        // value here; FlowModel will explicitly ignore invalid sections.
        const auto integrals = complete ? positiveGapIntegralsYZ(left, right, sections.gapSign)
                                        : GapIntegrals{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        sections.area[i] = integrals.area;
        sections.gapCubed[i] = integrals.cubed;
    }
}
