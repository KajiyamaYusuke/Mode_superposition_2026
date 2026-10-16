#include "SurfaceLoad.h"

void SurfaceLoad::resize(const Geometry& left, const Geometry& right) {
    fxL.assign(left.nsurfl, std::vector<double>(left.nsurfz, 0.0));
    fyL.assign(left.nsurfl, std::vector<double>(left.nsurfz, 0.0));
    fzL.assign(left.nsurfl, std::vector<double>(left.nsurfz, 0.0));
    fxR.assign(right.nsurfl, std::vector<double>(right.nsurfz, 0.0));
    fyR.assign(right.nsurfl, std::vector<double>(right.nsurfz, 0.0));
    fzR.assign(right.nsurfl, std::vector<double>(right.nsurfz, 0.0));
}
