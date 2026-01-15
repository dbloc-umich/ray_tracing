/**
 * Using Gauss-Legendre Quadrature to evaluate the integral of f(x) over a domain. 
**/

#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

#include "GaussianQuadrature.h"

class GaussLegendreBase{
    // Dummy class used to store the nodes and weights as static members
    public:
    virtual ~GaussLegendreBase() = default;
    inline static const std::vector<std::vector<double>> _nodes =
        {{}, // starting at n = 1
         {0.5773502691896257},
         {0.7745966692414834},
         {0.8611363115940526, 0.3399810435848563},
         {0.9061798459386640, 0.5384693101056831},
         {0.9324695142031521, 0.6612093864662645, 0.2386191860831969},
         {0.9491079123427585, 0.7415311855993945, 0.4058451513773972}
        };
    inline static const std::vector<std::vector<double>> _weights =
        {{2}, // starting at n = 1
         {1},
         {0.5555555555555556, 0.8888888888888889},
         {0.3478548451374539, 0.6521451548625461},
         {0.2369268850561891, 0.4786286704993665, 0.5688888888888889},
         {0.1713244923791704, 0.3607615730481386, 0.4679139345726910},
         {0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694}
        };
};

template<int N = 1, int M = 1>
class GaussLegendre : public GaussianQuadrature<N, M>, GaussLegendreBase{
    public:
    using typename GaussianQuadrature<N, M>::DomainType;
    using typename Quadrature<N, M>::RangeType;
    using typename Quadrature<N, M>::Function;
    using typename GaussianQuadrature<N, M>::IntegrationDomain;

    GaussLegendre(unsigned short n): GaussianQuadrature<N, M>(n), GaussLegendreBase()
    {
        if (n < 1 || n > 7) throw std::domain_error("Error: Can only support 1 through 7 integration nodes.");
    }

    RangeType integrate(const Function& f, const IntegrationDomain& D) const override{
        if (!this->validIntegrationDomain(D)) throw std::invalid_argument("ERROR: Invalid integration domain.");

        RangeType sum;
        if constexpr (M == 1) sum = 0.0;
        else if constexpr (M == Eigen::Dynamic) sum = RangeType::Zero(M,1);
        else M = RangeType::Zero();

        if constexpr (N == 1){
            auto nodes = getNodes(D.first, D.second);
            auto weights = getWeights(D.first, D.second);
            for (auto i = 0; i < this->_n; i++) sum += weights[i] * f(nodes[i]);
        } else{
            constexpr unsigned short dim = (N != Eigen::Dynamic) ? N : D.size();
            std::size_t numNodes = 1; // total number of nodes in all dimensions;
            for (auto i = 0; i < dim; i++) numNodes *= this->_n;

            std::vector<DomainType> nodes(numNodes);
            if (N == Eigen::Dynamic){
                for (auto node: nodes) node.resize(dim);
            }
            std::vector<double> weights(numNodes, 1.0);
            this->getNodesAndWeights(0, dim, 0, numNodes, nodes, weights, D);
            for (std::size_t i = 0; i < numNodes; i++) sum += weights[i] * f(nodes[i]);
        }
        return sum;
    }

    protected:
    Eigen::VectorXd getNodes(double l, double u) const noexcept override{
        Eigen::VectorXd nodes(this->_n);
        for (Eigen::Index i = 0; i < this->_n/2; i++){
            nodes[2*i] = _nodes[this->_n-1][i] * (u-l)/2 + (u+l)/2;
            nodes[2*i+1] = -_nodes[this->_n-1][i] * (u-l)/2 + (u+l)/2;
        }
        if (this->_n%2 == 1) nodes[this->_n-1] = (u+l)/2; // odd power, include a zero root
        return nodes;
    }

    Eigen::VectorXd getWeights(double l, double u) const noexcept override{
        Eigen::VectorXd weights(this->_n);
        for (Eigen::Index i = 0; i < this->_n/2; i+=2){
            weights[2*i] = _weights[this->_n-1][i] * (u-l)/2;
            weights[2*i+1] = weights[2*i];
        }
        if (this->_n%2 == 1) weights[this->_n-1] = _weights[this->_n-1].back() * (u-l)/2; // odd power, include the weight of the zero root
        return weights;
    }
};

#endif