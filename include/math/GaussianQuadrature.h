#ifndef GAUSSIAN_QUADRATURE_H
#define GAUSSIAN_QUADRATURE_H

#include "Quadrature.h"

template<int N, int M>
class GaussianQuadrature : public Quadrature<N, M>{
    public:
    using typename Quadrature<N, M>::DomainType;
    using typename Quadrature<N, M>::IntegrationDomain;

    GaussianQuadrature(unsigned short n): Quadrature<N, M>(), _n(n) {}
    virtual ~GaussianQuadrature() = default;

    protected:
    unsigned short _n; // number of nodes to evaluate
    virtual Eigen::VectorXd getNodes(double l, double u) const noexcept = 0;
    virtual Eigen::VectorXd getWeights(double l, double u) const noexcept = 0;

    void getNodesAndWeights(unsigned short d, const unsigned short dim, std::size_t start, std::size_t stop,
                            std::vector<DomainType>& allNodes, std::vector<double>& allWeights, const IntegrationDomain& D) const{
        // Recursively form a Cartesian product matrix of nodes and their corresponding weights

        if (d == dim) return;
        double l, u;
        if (d == 0){
            l = this->eval(D[0]);
            u = this->eval(D[1]);
        } else if (d == 1){
            double param = allNodes[0][0];
            l = this->eval(D[2], param);
            u = this->eval(D[3], param);
        } else{
            Eigen::VectorXd param = allNodes[0](Eigen::seq(0, d));
            l = this->eval(D[2*d], param);
            u = this->eval(D[2*d+1], param);
        }
        
        auto nodes = getNodes(l, u);
        auto weights = getWeights(l, u);

        std::size_t toFill = (stop-start)/_n;
        for (auto i = start; i < stop; i++){
            allNodes[i][d] = nodes[(i-start)/toFill];
            allWeights[i] *= weights[(i-start)/toFill];
        }

        for (auto i = 0; i < _n; i++) getNodesAndWeights(d+1, dim, start+i*toFill, start+(i+1)*toFill, allNodes, allWeights, D);
    }
};

#endif