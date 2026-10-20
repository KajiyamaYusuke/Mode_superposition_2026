#include "ChannelSections.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace {

// Millimetres. These tolerances merge shared cut vertices, not distinct branches.
constexpr double cutTolerance = 1.0e-9;
constexpr double minimumCoverage = 0.90;
constexpr int planeCount = 81;
struct CutPoint { double y, z; };
struct CutEdge { int a, b; };
struct Contour {
    std::vector<CutPoint> nodes;
    std::vector<CutEdge> edges;
    std::map<std::pair<long long, long long>, std::vector<int>> buckets;
    std::set<std::pair<int, int>> edgeKeys;
    bool coplanarFace = false;

    void clear() {
        nodes.clear(); edges.clear(); buckets.clear(); edgeKeys.clear();
        coplanarFace = false;
    }
    int node(CutPoint p) {
        const long long iy = std::llround(p.y / cutTolerance);
        const long long iz = std::llround(p.z / cutTolerance);
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                const auto it = buckets.find({iy + dy, iz + dz});
                if (it == buckets.end()) continue;
                for (int id : it->second) {
                    if (std::hypot(nodes[id].y - p.y, nodes[id].z - p.z)
                        <= cutTolerance) return id;
                }
            }
        }
        const int id = static_cast<int>(nodes.size());
        nodes.push_back(p);
        buckets[{iy, iz}].push_back(id);
        return id;
    }
    void add(CutPoint a, CutPoint b) {
        const int ia = node(a), ib = node(b);
        if (ia != ib && edgeKeys.insert(std::minmax(ia, ib)).second)
            edges.push_back({ia, ib});
    }
};

struct SurfaceCuts {
    const Geometry* geometry = nullptr;
    int rows = 0, columns = 0;
    std::vector<std::array<int, 3>> triangles;
    std::array<Contour, planeCount> cuts;
    double referenceZLow = 0.0, referenceZHigh = 0.0;

    void prepare(const Geometry& geom) {
        if (geometry == &geom && rows == geom.nsurfl && columns == geom.nsurfz) return;
        geometry = &geom; rows = geom.nsurfl; columns = geom.nsurfz;
        if (rows < 2 || columns < 2 || static_cast<int>(geom.surfp.size()) < rows)
            throw std::runtime_error("channel surface grid is incomplete");
        triangles.clear();
        referenceZLow = std::numeric_limits<double>::infinity();
        referenceZHigh = -referenceZLow;
        for (int i = 0; i < rows; ++i) {
            if (static_cast<int>(geom.surfp[i].size()) < columns)
                throw std::runtime_error("channel surface row is incomplete");
            for (int j = 0; j < columns; ++j) {
                const int id = geom.surfp[i][j];
                if (id < 0 || id >= static_cast<int>(geom.points.size()))
                    throw std::runtime_error("channel surface has an invalid node index");
                const double z = geom.points[id].z;
                if (!std::isfinite(z))
                    throw std::runtime_error("non-finite reference AP coordinate");
                referenceZLow = std::min(referenceZLow, z);
                referenceZHigh = std::max(referenceZHigh, z);
            }
        }
        // All nsurfl rows, not the nxsup pressure-load subset. Match the
        // diagnostic mesh diagonal p00--p11 without importing its solver.
        for (int i = 0; i + 1 < rows; ++i) {
            for (int j = 0; j + 1 < columns; ++j) {
                const int a = geom.surfp[i][j], b = geom.surfp[i + 1][j];
                const int c = geom.surfp[i + 1][j + 1], d = geom.surfp[i][j + 1];
                triangles.push_back({a, b, c});
                triangles.push_back({a, c, d});
            }
        }
    }

    void build(const Geometry& geom, const State& state, const std::vector<double>& planes) {
        prepare(geom);
        for (auto& cut : cuts) cut.clear();
        for (const auto& ids : triangles) {
            std::array<Point, 3> p;
            for (int v = 0; v < 3; ++v) {
                if (ids[v] >= static_cast<int>(state.disp.size()))
                    throw std::runtime_error("channel surface state is incomplete");
                const auto& u = state.disp[ids[v]];
                p[v] = {u.ux, u.uy, u.uz}; // Absolute deformed coordinates.
                if (!std::isfinite(u.ux) || !std::isfinite(u.uy) || !std::isfinite(u.uz))
                    throw std::runtime_error("non-finite channel surface node " + std::to_string(ids[v]));
                // Also bound spatial-hash conversion well inside long long.
                if (std::abs(u.uy) > 1.0e8 || std::abs(u.uz) > 1.0e8)
                    throw std::runtime_error("channel surface coordinates exceed cut-index range");
            }
            const double low = std::min({p[0].x, p[1].x, p[2].x});
            const double high = std::max({p[0].x, p[1].x, p[2].x});
            // Bucket by each triangle's x extent: never scan all triangles at
            // every (x,z) sample.
            const auto first = std::lower_bound(planes.begin(), planes.end(), low - cutTolerance);
            const auto last = std::upper_bound(planes.begin(), planes.end(), high + cutTolerance);
            for (auto it = first; it != last; ++it) {
                auto& cut = cuts[it - planes.begin()];
                std::vector<CutPoint> hit;
                auto addHit = [&](CutPoint q) {
                    for (const auto& h : hit)
                        if (std::hypot(h.y - q.y, h.z - q.z) <= cutTolerance) return;
                    hit.push_back(q);
                };
                int onPlane = 0;
                for (int v = 0; v < 3; ++v) {
                    if (std::abs(p[v].x - *it) <= cutTolerance) {
                        ++onPlane;
                        addHit({p[v].y, p[v].z});
                    }
                }
                if (onPlane == 3) {
                    // A two-dimensional coplanar face is not a unique
                    // one-dimensional contour. Report it if this plane is used.
                    cut.coplanarFace = true;
                    for (int v = 0; v < 3; ++v)
                        cut.add({p[v].y, p[v].z}, {p[(v + 1) % 3].y, p[(v + 1) % 3].z});
                    continue;
                }
                for (int v = 0; v < 3; ++v) {
                    const auto& a = p[v]; const auto& b = p[(v + 1) % 3];
                    const double da = a.x - *it, db = b.x - *it;
                    if ((da < -cutTolerance && db > cutTolerance)
                        || (db < -cutTolerance && da > cutTolerance)) {
                        const double s = (*it - a.x) / (b.x - a.x);
                        addHit({a.y + s * (b.y - a.y), a.z + s * (b.z - a.z)});
                    }
                }
                if (hit.size() == 2) cut.add(hit[0], hit[1]);
            }
        }
    }
};

double cross(CutPoint a, CutPoint b, CutPoint c) {
    return (b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y);
}

void validateContour(const Contour& contour) {
    if (contour.coplanarFace) throw std::runtime_error("coplanar surface face: contour is not unique");
    if (contour.edges.empty()) throw std::runtime_error("missing surface contour");
    std::vector<std::vector<int>> adjacent(contour.nodes.size());
    for (const auto& edge : contour.edges) {
        adjacent[edge.a].push_back(edge.b); adjacent[edge.b].push_back(edge.a);
    }
    for (const auto& neighbours : adjacent)
        if (neighbours.size() > 2) throw std::runtime_error("branched surface contour");
    std::vector<bool> visited(contour.nodes.size(), false);
    std::vector<int> stack{contour.edges.front().a};
    visited[stack.front()] = true;
    while (!stack.empty()) {
        const int a = stack.back(); stack.pop_back();
        for (int b : adjacent[a]) if (!visited[b]) { visited[b] = true; stack.push_back(b); }
    }
    for (size_t i = 0; i < adjacent.size(); ++i)
        if (!adjacent[i].empty() && !visited[i])
            throw std::runtime_error("multiple independent surface contours");

    // Bounding-box rejection keeps the self-crossing check local to
    // overlapping yz segments. Shared endpoints alone are not crossings.
    for (size_t i = 0; i < contour.edges.size(); ++i) {
        const auto e = contour.edges[i];
        const auto a = contour.nodes[e.a], b = contour.nodes[e.b];
        for (size_t j = i + 1; j < contour.edges.size(); ++j) {
            const auto f = contour.edges[j];
            const auto c = contour.nodes[f.a], d = contour.nodes[f.b];
            if (std::max(a.y, b.y) + cutTolerance < std::min(c.y, d.y)
                || std::max(c.y, d.y) + cutTolerance < std::min(a.y, b.y)
                || std::max(a.z, b.z) + cutTolerance < std::min(c.z, d.z)
                || std::max(c.z, d.z) + cutTolerance < std::min(a.z, b.z)) continue;
            const double eps = cutTolerance * std::max({1.0,
                std::hypot(b.y - a.y, b.z - a.z), std::hypot(d.y - c.y, d.z - c.z)});
            const double ac = cross(a, b, c), ad = cross(a, b, d);
            const double ca = cross(c, d, a), cb = cross(c, d, b);
            const bool shared = e.a == f.a || e.a == f.b || e.b == f.a || e.b == f.b;
            if (shared) {
                if (std::abs(ac) <= eps && std::abs(ad) <= eps) {
                    const bool useY = std::abs(b.y - a.y) >= std::abs(b.z - a.z);
                    const double a0 = useY ? a.y : a.z, a1 = useY ? b.y : b.z;
                    const double b0 = useY ? c.y : c.z, b1 = useY ? d.y : d.z;
                    if (std::min(std::max(a0, a1), std::max(b0, b1))
                        - std::max(std::min(a0, a1), std::min(b0, b1)) > cutTolerance)
                        throw std::runtime_error("self-overlapping surface contour");
                }
                continue;
            }
            auto opposite = [eps](double u, double v) {
                return (u <= eps && v >= -eps) || (v <= eps && u >= -eps);
            };
            if (opposite(ac, ad) && opposite(ca, cb))
                throw std::runtime_error("self-intersecting surface contour");
        }
    }
}

bool hasAtZ(const Contour& contour, double z) {
    for (const auto& edge : contour.edges) {
        const auto a = contour.nodes[edge.a], b = contour.nodes[edge.b];
        if (z > std::min(a.z, b.z) && z < std::max(a.z, b.z)) return true;
    }
    return false;
}

double yOnEdge(const Contour& contour, int edgeId, double z) {
    const auto edge = contour.edges[edgeId];
    const auto a = contour.nodes[edge.a], b = contour.nodes[edge.b];
    return a.y + (z - a.z) / (b.z - a.z) * (b.y - a.y);
}

// At a knot include horizontal edges and deduplicate real shared vertices.
double envelopeAtKnot(const Contour& contour, double z, bool takeMaximum) {
    std::vector<double> ys;
    for (const auto& edge : contour.edges) {
        const auto a = contour.nodes[edge.a], b = contour.nodes[edge.b];
        if (z < std::min(a.z, b.z) - cutTolerance
            || z > std::max(a.z, b.z) + cutTolerance) continue;
        if (std::abs(b.z - a.z) <= cutTolerance) {
            ys.push_back(a.y); ys.push_back(b.y);
        } else {
            const double s = std::clamp((z - a.z) / (b.z - a.z), 0.0, 1.0);
            ys.push_back(a.y + s * (b.y - a.y));
        }
    }
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end(), [](double a, double b) {
        return std::abs(a - b) <= cutTolerance;
    }), ys.end());
    if (ys.size() > 2) throw std::runtime_error("three or more distinct intersections at z=" + std::to_string(z));
    return ys.empty() ? std::numeric_limits<double>::quiet_NaN()
                      : (takeMaximum ? ys.back() : ys.front());
}

int envelopeEdge(const Contour& contour, double z, bool takeMaximum) {
    int selected = -1;
    std::vector<double> ys;
    for (int i = 0; i < static_cast<int>(contour.edges.size()); ++i) {
        const auto edge = contour.edges[i];
        const auto a = contour.nodes[edge.a], b = contour.nodes[edge.b];
        if (!(z > std::min(a.z, b.z) && z < std::max(a.z, b.z))) continue;
        const double y = yOnEdge(contour, i, z);
        bool duplicate = false;
        for (double other : ys) if (std::abs(y - other) <= cutTolerance) duplicate = true;
        if (!duplicate) ys.push_back(y);
        if (selected < 0 || (takeMaximum ? y > yOnEdge(contour, selected, z)
                                        : y < yOnEdge(contour, selected, z))) selected = i;
    }
    if (ys.size() > 2) throw std::runtime_error("three or more distinct intersections at z=" + std::to_string(z));
    return selected;
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

} // namespace

void ChannelSections::resize(const Geometry& left, const Geometry& right, int nSections) {
    if (nSections != planeCount) throw std::runtime_error("fixed-x channel requires 81 sections");
    // Legacy material sigma/anchor fields are deliberately inactive.
    sigma.clear();
    referenceConstrictionI = -1;
    referenceConstrictionSigma = std::numeric_limits<double>::quiet_NaN();
    x.resize(planeCount);
    for (int i = 0; i < planeCount; ++i)
        x[i] = i <= 20 ? 0.5 * i : 10.0 + 0.05 * (i - 20);
    area.assign(planeCount, std::numeric_limits<double>::quiet_NaN());
    gapCubed = area;
    valid.assign(planeCount, false);
    double initialGap = 0.0;
    int count = 0;
    for (int i = 0; i < std::min(left.nsurfl, right.nsurfl); ++i) {
        for (int j = 0; j < std::min(left.nsurfz, right.nsurfz); ++j) {
            const int il = left.surfp.at(i).at(j), ir = right.surfp.at(i).at(j);
            if (il < 0 || ir < 0) continue;
            initialGap += right.points.at(ir).y - left.points.at(il).y;
            ++count;
        }
    }
    if (count == 0 || !std::isfinite(initialGap) || initialGap == 0.0)
        throw std::runtime_error("initial left/right gap direction is undefined");
    gapSign = initialGap < 0.0 ? -1.0 : 1.0;
}

void ChannelSectionBuilder::build(const Geometry& leftGeom, const State& leftState,
                                  const Geometry& rightGeom, const State& rightState,
                                  ChannelSections& sections) const {
    if (sections.x.size() != planeCount) throw std::runtime_error("fixed-x section grid is uninitialized");
    // The application builds one channel serially. Cache only static mesh
    // connectivity and scratch capacity; no physical state is cached.
    static thread_local SurfaceCuts left, right;
    static thread_local const ChannelSections* loggedSections = nullptr;
    static thread_local int loggedLast = -2;
    left.build(leftGeom, leftState, sections.x);
    right.build(rightGeom, rightState, sections.x);
    const double refLow = std::max(left.referenceZLow, right.referenceZLow);
    const double refHigh = std::min(left.referenceZHigh, right.referenceZHigh);
    if (!(refHigh > refLow)) throw std::runtime_error("reference AP span is empty");
    std::fill(sections.valid.begin(), sections.valid.end(), false);
    std::fill(sections.area.begin(), sections.area.end(), std::numeric_limits<double>::quiet_NaN());
    std::fill(sections.gapCubed.begin(), sections.gapCubed.end(), std::numeric_limits<double>::quiet_NaN());

    int lastValid = -1;
    double finalCoverage = 0.0;
    for (int k = 0; k < planeCount; ++k) {
        try {
            const auto& l = left.cuts[k]; const auto& r = right.cuts[k];
            std::vector<double> knots;
            for (const auto& p : l.nodes) knots.push_back(p.z);
            for (const auto& p : r.nodes) knots.push_back(p.z);
            for (int j = 0; refLow + 0.1 * j < refHigh; ++j) knots.push_back(refLow + 0.1 * j);
            knots.push_back(refHigh);
            std::sort(knots.begin(), knots.end());
            knots.erase(std::unique(knots.begin(), knots.end(), [](double a, double b) {
                return std::abs(a - b) <= cutTolerance;
            }), knots.end());
            double coveredLength = 0.0;
            for (size_t j = 1; j < knots.size(); ++j) {
                const double mid = 0.5 * (knots[j - 1] + knots[j]);
                if (hasAtZ(l, mid) && hasAtZ(r, mid)) coveredLength += knots[j] - knots[j - 1];
            }
            finalCoverage = coveredLength / (refHigh - refLow);
            if (finalCoverage < minimumCoverage) break;

            validateContour(l); validateContour(r);
            const bool maxLeft = sections.gapSign > 0.0;
            std::vector<double> leftAtKnots, rightAtKnots;
            leftAtKnots.reserve(knots.size()); rightAtKnots.reserve(knots.size());
            for (double z : knots) {
                leftAtKnots.push_back(envelopeAtKnot(l, z, maxLeft));
                rightAtKnots.push_back(envelopeAtKnot(r, z, !maxLeft));
            }
            double area = 0.0, cubed = 0.0;
            for (size_t j = 1; j < knots.size(); ++j) {
                const double z0 = knots[j - 1], z1 = knots[j], mid = 0.5 * (z0 + z1);
                const int il = envelopeEdge(l, mid, maxLeft), ir = envelopeEdge(r, mid, !maxLeft);
                const double l0 = il >= 0 ? yOnEdge(l, il, z0) : 0.0;
                const double l1 = il >= 0 ? yOnEdge(l, il, z1) : 0.0;
                const double r0 = ir >= 0 ? yOnEdge(r, ir, z0) : 0.0;
                const double r1 = ir >= 0 ? yOnEdge(r, ir, z1) : 0.0;
                if ((il >= 0 && (std::abs(l0 - leftAtKnots[j - 1]) > cutTolerance
                                 || std::abs(l1 - leftAtKnots[j]) > cutTolerance))
                    || (ir >= 0 && (std::abs(r0 - rightAtKnots[j - 1]) > cutTolerance
                                    || std::abs(r1 - rightAtKnots[j]) > cutTolerance)))
                    throw std::runtime_error("discontinuous inner envelope / branch switch");
                if (il < 0 || ir < 0) continue; // Never bridge absent AP intervals.
                const double h0 = sections.gapSign * (r0 - l0);
                const double h1 = sections.gapSign * (r1 - l1);
                positiveLinearIntegrals(h0, h1, z1 - z0, area, cubed);
            }
            if (!std::isfinite(area) || !std::isfinite(cubed))
                throw std::runtime_error("non-finite section integrals");
            // Observed inner-envelope area, not a proven full airway section.
            sections.area[k] = area;
            sections.gapCubed[k] = cubed;
            sections.valid[k] = true;
            lastValid = k;
        } catch (const std::exception& error) {
            throw std::runtime_error("section=" + std::to_string(k)
                + " x_mm=" + std::to_string(sections.x[k]) + ": " + error.what());
        }
    }
    if (loggedSections != &sections || loggedLast != lastValid) {
        std::cout << "[ChannelSections] integration_end_x_mm="
                  << (lastValid >= 0 ? sections.x[lastValid] : std::numeric_limits<double>::quiet_NaN());
        if (lastValid + 1 < planeCount)
            std::cout << " first_invalid_x_mm=" << sections.x[lastValid + 1]
                      << " AP_coverage=" << finalCoverage << " (required >=0.90)";
        std::cout << "; contiguous inlet prefix only, no downstream restart" << std::endl;
        loggedSections = &sections; loggedLast = lastValid;
    }
    if (lastValid < 1) throw std::runtime_error("section=0: insufficient AP coverage for an inlet-connected channel");
}
