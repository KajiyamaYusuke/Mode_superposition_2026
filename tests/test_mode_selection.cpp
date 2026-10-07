#include "ModeData.h"

#include <cassert>

int main() {
    Geometry geometry;
    geometry.nPoints = 2;

    ModeData defaults;
    defaults.initialize(5, geometry);
    assert(defaults.nModes == 5);
    for (int m = 0; m < 5; ++m) assert(defaults.sourceModeIndex(m) == m);

    ModeData fifthOnly;
    fifthOnly.initialize(10, geometry, {5});
    assert(fifthOnly.nModes == 1);
    assert(fifthOnly.sourceModeIndex(0) == 4);
    assert(fifthOnly.modes.size() == 1);

    ModeData selected;
    selected.initialize(10, geometry, {1, 5, 8});
    assert(selected.nModes == 3);
    assert(selected.sourceModeIndex(0) == 0);
    assert(selected.sourceModeIndex(1) == 4);
    assert(selected.sourceModeIndex(2) == 7);
}
