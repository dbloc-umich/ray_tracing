#include "ExplicitMidpoint.h"

template<int M>
ODEStatus ExplicitMidpoint<M>::solve(const Function& f, const RangeType& ic) const noexcept{
    // Calculates the solution array size
    Eigen::Index steps = this->_t.size();
    if (steps == 0) return ODEStatus::InvalidArgument;
    Eigen::Index m;
    if constexpr(M == Eigen::Dynamic) m = ic.size();
    else m = M;

    if (this->_u.rows() != m || this->_u.cols() != steps) return ODEStatus::InvalidArgument;
    this->_u.col(0) = ic; // Initialize the "state" vector

    // Solve
    for (Eigen::Index n = 1; n < steps; n++){
        double dt = this->_t[n] - this->_t[n-1];
        RangeType k1, k2;
        if constexpr(M == 1){
            k1 = f(this->_t[n-1], this->_u[n-1]);
            k2 = f(this->_t[n-1] + 0.5*dt, this->_u[n-1] + 0.5*dt*k1);
            this->_u[n] = this->_u[n-1] + dt*k2;
        } else{
            k1 = f(this->_t[n-1], this->_u.col(n-1));
            k2 = f(this->_t[n-1] + 0.5*dt, this->_u.col(n-1) + 0.5*dt*k1);
            this->_u.col(n) = this->_u.col(n-1) + dt*k2;
        }

        if (!Eigen::isfinite(this->_u.col(n)).all()){
            this->_t = this->_t.head(n+1);
            this->_u = this->_u.topLeftCorner(m,n+1);
            return ODEStatus::FailureToSolve;
        }
    }
    return ODEStatus::Success;
}