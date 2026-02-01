#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include "ODE_IVPSolver.h"
#include "NewtonSolver.h"

template<int M>
class CrankNicolson: public ODE_IVPSolver<M>{
    public:
    using typename ODESolver<M>::RangeType;
    using typename ODESolver<M>::SolutionType;
    using typename ODE_IVPSolver<M>::Function;
    using SolverPointer = std::unique_ptr<NonlinearSolver<M>>;
    
    CrankNicolson(double t0, double t1, double dt, SolverPointer solver=std::make_unique<NewtonSolver<M>>()):
        ODE_IVPSolver<M>(t0, t1, dt),
        _nlSolver(std::move(solver))
    {}
    CrankNicolson(const Eigen::ArrayXd& t, SolverPointer solver=std::make_unique<NewtonSolver<M>>()):
        ODE_IVPSolver<M>(t),
        _nlSolver(std::move(solver))
    {}    

    ODEStatus solve(const Function& f, const RangeType& ic) const noexcept override{
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
                k2 = k1;
                auto k2Func = [&, this](const RangeType& k2){ return k2 - f(this->_t[n], this->_u[n-1] + 0.5*dt*(k1+k2)); };
                _nlSolver->setFunction(k2Func);
                auto status = _nlSolver->solve(k2);
                switch (status){
                    case (NLStatus::Success):
                        this->_u[n] = this->_u[n-1] + 0.5*dt*(k1+k2);
                        break;
                    default:
                        this->_t = this->_t.head(n+1);
                        this->_u = this->_u.topLeftCorner(m,n+1);
                        return ODEStatus::FailureToSolve;
                }
                
            } else{
                k1 = f(this->_t[n-1], this->_u.col(n-1));
                k2 = k1;
                auto k2Func = [&, this](const RangeType& k2){ return k2 - f(this->_t[n], this->_u.col(n-1) + 0.5*dt*(k1+k2)); };
                _nlSolver->setFunction(k2Func);
                auto status = _nlSolver->solve(k2);
                switch (status){
                    case (NLStatus::Success):
                        this->_u.col(n) = this->_u.col(n-1) + 0.5*dt*(k1+k2);
                        break;
                    default:
                        this->_t = this->_t.head(n+1);
                        this->_u = this->_u.topLeftCorner(m,n+1);
                        return ODEStatus::FailureToSolve;
                }
            }

            if (!Eigen::isfinite(this->_u.col(n)).all()){
                this->_t = this->_t.head(n+1);
                this->_u = this->_u.topLeftCorner(m,n+1);
                return ODEStatus::FailureToSolve;
            }
        }
        return ODEStatus::Success;
    }

    protected:
    SolverPointer _nlSolver;   
};

#endif