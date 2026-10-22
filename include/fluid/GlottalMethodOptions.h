#pragma once
#include <map>
#include <stdexcept>
#include <string>

struct GlottalMethodOptions {
    std::string profile = "modern";
    bool area = false;
    bool flowCoefficients = false;
    bool wallViscosity = false;
    bool separation = false;
    bool postSeparation = false;
    bool surfacePressureLoad = false;

    bool any() const {
        return area || flowCoefficients || wallViscosity || separation
            || postSeparation || surfacePressureLoad;
    }
    void configure(const std::string& name, const std::map<std::string, bool>& overrides) {
        if (name != "modern" && name != "legacy_fluid_all" && name != "custom")
            throw std::runtime_error("unknown glottalMethodProfile: " + name);
        *this = GlottalMethodOptions{};
        profile = name;
        const bool old = name == "legacy_fluid_all";
        area = flowCoefficients = wallViscosity = separation = postSeparation = surfacePressureLoad = old;
        std::map<std::string, bool*> flags{{"legacyarea", &area},
            {"legacyflowcoefficients", &flowCoefficients}, {"legacywallviscosity", &wallViscosity},
            {"legacyseparation", &separation}, {"legacypostseparation", &postSeparation},
            {"legacysurfacepressureload", &surfacePressureLoad}};
        for (const auto& entry : overrides) {
            if (entry.first == "legacycontact" || entry.first == "legacyintegrator") {
                if (entry.second) throw std::runtime_error(entry.first + " is not implemented; keep phase-1 contact/Newmark");
                continue;
            }
            if (!flags.count(entry.first)) throw std::runtime_error("unknown legacy flag: " + entry.first);
            if ((name == "modern" && entry.second) || (old && !entry.second))
                throw std::runtime_error("profile/flag conflict; use glottalMethodProfile=custom");
            *flags.at(entry.first) = entry.second;
        }
        if (surfacePressureLoad && !area)
            throw std::runtime_error("legacy material-row pressure load requires legacyArea; change area+mapping together");
    }
};
