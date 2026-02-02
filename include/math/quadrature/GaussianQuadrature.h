#ifndef GAUSSIAN_QUADRATURE_H
#define GAUSSIAN_QUADRATURE_H

#include "Quadrature.h"

template<int N = 1, int M = 1>
class GaussianQuadrature : public Quadrature<N, M>{
    public:
    using typename Quadrature<N, M>::DomainType;
    using typename Quadrature<N, M>::RangeType;
    using typename Quadrature<N, M>::Function;
    using typename Quadrature<N, M>::IntegrationDomain;

    GaussianQuadrature(std::size_t n): Quadrature<N, M>(), _n(n) {}
    virtual ~GaussianQuadrature() = default;

    virtual RangeType integrate(const Function& f, const IntegrationDomain& D) const override;

    protected:
    std::size_t _n; // number of nodes to evaluate
    virtual std::size_t n() const noexcept = 0; // number of actual nodes per dimension
    virtual Eigen::VectorXd getNodes(double l, double u) const = 0;
    virtual Eigen::VectorXd getWeights(double l, double u) const = 0;
    void getNodesAndWeights(std::size_t d, const std::size_t dim, std::size_t start, std::size_t stop,
                            const IntegrationDomain& D, std::vector<DomainType>& allNodes, std::vector<double>& allWeights) const;
};

#include "GaussianQuadrature.tpp"
#endif