#pragma once

// All closure decisions use this single SI threshold.  Area in the existing
// surface/section code is expressed in mm^2, hence the paired conversion.
constexpr double kAreaClosedM2 = 1.0e-8;
constexpr double kAreaClosedMm2 = kAreaClosedM2 * 1.0e6;
