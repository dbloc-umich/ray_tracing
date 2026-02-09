#include "ODESolver.h"

template<int M>
void ODESolver<M>::makeGrid(double t0, double t1, double dt) const noexcept{
    assert(dt != 0 && t1 != t0);
    assert((t1 > t0 && dt > 0) || (t1 < t0 && dt < 0));
    Eigen::Index steps = ceil((t1-t0)/dt);
    _t.resize(steps+1);
    _t[0] = t0;
    _t[steps] = t1;
    for (Eigen::Index i = 1; i < steps; i++) _t[i] = _t[i-1] + dt;
    _u.resize(M == Eigen::Dynamic ? 1 : M, steps+1);
}

template<int M>
void ODESolver<M>::makeGrid(const Eigen::ArrayXd& t) const noexcept{
    assert(std::is_sorted(t.cbegin(), t.cend(), [](double a, double b){ return a < b; } ));
    _t = t;
    _u.resize(M == Eigen::Dynamic ? 1 : M, _t.size());
}