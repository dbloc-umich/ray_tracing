#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "Eigen/Dense"
#include <functional>
#include <variant>
#include <vector>

template<int N = 1, int M = 1>
class Quadrature{
    static_assert(N == Eigen::Dynamic || N >= 1);
    static_assert(M == Eigen::Dynamic || M >= 1);
    public:
    using DomainType = std::conditional_t<N == 1, double, Eigen::Matrix<double, N, 1>>;
    using RangeType = std::conditional_t<M == 1, double, Eigen::Matrix<double, M, 1>>;
    using Function = std::function<RangeType(const DomainType&)>;

    using IFunc1D = std::function<double(double)>; // 1D function as an integration bound
    using IFuncMD = std::function<double(const Eigen::VectorXd&)>; // MD function as an integration bound
    using IntegrationBound = std::variant<double, IFunc1D, IFuncMD>;
    using IntegrationDomain = std::vector<IntegrationBound>;

    virtual ~Quadrature() = default;
    virtual RangeType integrate(const Function&, const IntegrationDomain&) const = 0;

    protected:
    bool validIntegrationDomain(const IntegrationDomain& D) const noexcept;
    double eval(const IntegrationBound& B) const;
    double eval(const IntegrationBound& B, double param) const;
    double eval(const IntegrationBound& B, const Eigen::VectorXd& param) const;
};

#include "Quadrature.tpp"
#endif