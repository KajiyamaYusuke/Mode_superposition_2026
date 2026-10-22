#include "LegacyGlottalMethods.h"
#include <cmath>
#include <stdexcept>

LegacyGlottal::Coefficients LegacyGlottal::coefficients(double rho, double mu, double beta,
                                                       double lis, double lg, double area) {
    if (!(rho>0 && mu>=0 && beta>0 && lis>0 && lg>0 && area>0)
        || !std::isfinite(rho+mu+beta+lis+lg+area))
        throw std::runtime_error("invalid legacy flow coefficient input (all lengths/areas must be SI)");
    return {rho*.5*lis/area, 12.*mu*lg*lg*lis/(area*area*area), beta*rho/(area*area)};
}
