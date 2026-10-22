// Geometry preparation only. Reuse unchanged production initialization,
// mode normalization/mirroring and section rules; never call Simulation::run.
#include "../../src/fluid/channelSections.cpp"
#include "Simulation.h"
#include "SeparationSelection.h"
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>

template<class T> void put(std::ofstream& file, T value) {
    file.write(reinterpret_cast<const char*>(&value), sizeof(value));
    if (!file) throw std::runtime_error("binary preparation write failed");
}

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: basis params native.bin new-preparation-directory");
        const fs::path target = argv[3];
        if (fs::exists(target / "basis.json")) throw std::runtime_error("refusing to overwrite basis");
        fs::create_directories(target);
        Simulation sim;
        sim.initialize(argv[1]);
        if (sim.mdataL.nModes != 20 || sim.mdataR.nModes != 20)
            throw std::runtime_error("expected all 20 modes on both folds");
        std::ofstream reference(target / "reference.bin", std::ios::binary);
        std::ofstream modes(target / "modes.bin", std::ios::binary);
        std::ofstream connectivity(target / "connectivity.bin", std::ios::binary);
        std::ofstream offsets(target / "offsets.bin", std::ios::binary);
        std::ofstream types(target / "types.bin", std::ios::binary);
        std::int64_t pointBase = 0, connectionEnd = 0;
        for (int side = 0; side < 2; ++side) {
            const auto& geom = side == 0 ? sim.geomL : sim.geomR;
            const auto& md = side == 0 ? sim.mdataL : sim.mdataR;
            for (const auto& point : geom.points)
                for (double value : {point.x, point.y, point.z}) put(reference, value);
            for (int m = 0; m < md.nModes; ++m)
                for (const auto& point : md.modes[m])
                    for (double value : {point.ux, point.uy, point.uz}) put(modes, value);
            if (geom.types.size() != geom.connect.size()) throw std::runtime_error("cell type count mismatch");
            for (std::size_t c = 0; c < geom.connect.size(); ++c) {
                for (int id : geom.connect[c]) {
                    if (id < 0 || id >= geom.nPoints) throw std::runtime_error("invalid original connectivity");
                    put(connectivity, pointBase + id);
                    ++connectionEnd;
                }
                put(offsets, connectionEnd);
                put(types, static_cast<std::uint8_t>(geom.types[c]));
            }
            pointBase += geom.nPoints;
        }
        std::ofstream basis(target / "basis.json");
        basis << std::setprecision(17)
              << "{\n\"points_per_fold\": [" << sim.geomL.nPoints << ',' << sim.geomR.nPoints
              << "],\n\"cells_per_fold\": [" << sim.geomL.nCells << ',' << sim.geomR.nCells
              << "],\n\"modes_per_fold\": 20,\n\"dt_s\": " << sim.params.dt
              << ",\n\"mass_kg\": " << sim.params.mass
              << ",\n\"reference_order\": \"fold L/R, point, xyz; little endian float64 mm\","
              << "\n\"mode_order\": \"fold L/R, mode, point, xyz; little endian float64\","
              << "\n\"connectivity\": \"int64, left cells then right cells with point offset\","
              << "\n\"normalization\": \"unchanged Simulation::initialize and ModeData::normalizeModes\","
              << "\n\"coupled_steps_executed\": 0\n}\n";
        if (!basis) throw std::runtime_error("basis metadata write failed");
        std::ifstream native(argv[2], std::ios::binary);
        std::ofstream checks(target / "geometry_checks.csv");
        checks << "step,time_s,min_index,x_min_mm,Amin_mm2\n" << std::setprecision(17);
        // Only five control states, not a time integration or a full-dt mesh export.
        for (int step : {12000, 12480, 12792, 15001, 17000}) {
            native.seekg(std::streamoff(step) * 2 * 20 * 7 * 8);
            std::ofstream control(target / ("control_" + std::to_string(step) + ".bin"), std::ios::binary);
            for (int side = 0; side < 2; ++side) {
                auto& state = side == 0 ? sim.stateL : sim.stateR;
                const auto& geom = side == 0 ? sim.geomL : sim.geomR;
                const auto& md = side == 0 ? sim.mdataL : sim.mdataR;
                state.qf.resize(20);
                state.qfdot.resize(20);
                for (int m = 0; m < 20; ++m) {
                    std::array<double, 7> row{};
                    native.read(reinterpret_cast<char*>(row.data()), sizeof(row));
                    if (!native || !std::isfinite(row[4] + row[5])) throw std::runtime_error("invalid native history");
                    state.qf[m] = row[4];
                    state.qfdot[m] = row[5];
                }
                state.mode2uf(geom, md, step);
                state.disp = state.predictedDisp;
                for (const auto& point : state.disp)
                    for (double value : {point.ux, point.uy, point.uz}) put(control, value);
            }
            ChannelSections sections;
            sections.resize(sim.geomL, sim.geomR, 50);
            ChannelSectionBuilder().build(sim.geomL, sim.stateL, sim.geomR, sim.stateR, sections);
            const auto chosen = selectSeparationByAreaRatio(sections.area, sections.valid, 1., .01);
            checks << step << ',' << step * sim.params.dt << ',' << chosen.minIndex << ','
                   << sections.x[chosen.minIndex] << ',' << chosen.minAreaMm2 << '\n';
        }
        std::cout << "PASS geometry basis and five controls; no coupled steps executed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
