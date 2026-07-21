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

double polygonAreaYZ(const std::vector<SectionPoint>& left,
                     const std::vector<SectionPoint>& right) {
    std::vector<SectionPoint> polygon;
    polygon.reserve(left.size() + right.size());
    polygon.insert(polygon.end(), left.begin(), left.end());
    for (auto it = right.rbegin(); it != right.rend(); ++it) polygon.push_back(*it);
    if (polygon.size() < 3) return 0.0;

    double twiceArea = 0.0;
    for (size_t k = 0; k < polygon.size(); ++k) {
        const auto& a = polygon[k];
        const auto& b = polygon[(k + 1) % polygon.size()];
        twiceArea += a.y * b.z - b.y * a.z;
    }
    return 0.5 * std::abs(twiceArea);
}

void buildSurfaceAngles(const Geometry& geom, const State& state,
                        std::vector<std::vector<std::vector<double>>>& degree) {
    for (auto& component : degree)
        for (auto& row : component) std::fill(row.begin(), row.end(), 0.0);

    for (int i = 1; i < geom.nxsup - 1; ++i) {
        for (int j = 1; j < geom.nsurfz - 1; ++j) {
            const int il = geom.surfp[i - 1][j], ir = geom.surfp[i + 1][j];
            const int jd = geom.surfp[i][j - 1], ju = geom.surfp[i][j + 1];
            if (il < 0 || ir < 0 || jd < 0 || ju < 0) continue;
            const double dx = 0.5 * (state.disp[ir].ux - state.disp[il].ux);
            const double dyx = 0.5 * (state.disp[ir].uy - state.disp[il].uy);
            const double dz = 0.5 * (state.disp[ju].uz - state.disp[jd].uz);
            const double dyz = 0.5 * (state.disp[ju].uy - state.disp[jd].uy);
            if (std::abs(dx) > 1e-12) degree[0][i][j] = std::atan(dyx / dx);
            if (std::abs(dz) > 1e-12) degree[1][i][j] = std::atan(dyz / dz);
        }
    }
}
} // namespace

void ChannelSections::resize(const Geometry& left, const Geometry& right, int nSections) {
    const int nx = nSections;
    const int nzL = left.nsurfz;
    const int nzR = right.nsurfz;
    sigma.assign(nx, 0.0);
    x.assign(nx, 0.0);
    area.assign(nx, 0.0);
    valid.assign(nx, false);
    for (int k = 0; k < nx; ++k) sigma[k] = nx > 1 ? static_cast<double>(k) / (nx - 1) : 0.0;
    degreeL.assign(2, std::vector<std::vector<double>>(left.nxsup, std::vector<double>(nzL, 0.0)));
    degreeR.assign(2, std::vector<std::vector<double>>(right.nxsup, std::vector<double>(nzR, 0.0)));
}

void ChannelSectionBuilder::build(const Geometry& leftGeom, const State& leftState,
                                  const Geometry& rightGeom, const State& rightState,
                                  ChannelSections& sections) const {
    const int nx = static_cast<int>(sections.sigma.size());
    const int nz = std::min(leftGeom.nsurfz, rightGeom.nsurfz);
    if (nx < 2 || nz < 2) return;

    // Determine a deformed material centreline from the central z line.  The
    // sigma labels stay fixed, while their physical x positions can stretch.
    const int jMid = nz / 2;
    auto centerXAt = [&](double sigma) {
        const double u = sigma * (std::min(leftGeom.nxsup, rightGeom.nxsup) - 1);
        const int i0 = std::max(0, std::min(static_cast<int>(u), std::min(leftGeom.nxsup, rightGeom.nxsup) - 2));
        const double a = u - i0;
        const int l0 = leftGeom.surfp[i0][jMid], l1 = leftGeom.surfp[i0 + 1][jMid];
        const int r0 = rightGeom.surfp[i0][jMid], r1 = rightGeom.surfp[i0 + 1][jMid];
        return 0.5 * ((1.0-a)*leftState.disp[l0].ux + a*leftState.disp[l1].ux
                    + (1.0-a)*rightState.disp[r0].ux + a*rightState.disp[r1].ux);
    };

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
            sections.x[i] = centerXAt(sections.sigma[i]);
            sections.valid[i] = false;
            sections.area[i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        sections.x[i] = std::clamp(centerXAt(sections.sigma[i]), commonLow, commonHigh);
        std::vector<SectionPoint> left, right;
        left.reserve(nz); right.reserve(nz);
        bool complete = true;
        for (int j = 0; j < nz; ++j) {
            SectionPoint pl{}, pr{};
            complete = intersectAtX(leftGeom, leftState, j, sections.x[i], i, pl)
                    && intersectAtX(rightGeom, rightState, j, sections.x[i], i, pr);
            if (!complete) break;
            left.push_back(pl); right.push_back(pr);
        }
        sections.valid[i] = complete;
        // Invalid is distinct from a genuine closed section.  Keep a benign
        // value here; FlowModel will explicitly ignore invalid sections.
        sections.area[i] = complete ? polygonAreaYZ(left, right)
                                    : std::numeric_limits<double>::quiet_NaN();
    }

    buildSurfaceAngles(leftGeom, leftState, sections.degreeL);
    buildSurfaceAngles(rightGeom, rightState, sections.degreeR);
}
