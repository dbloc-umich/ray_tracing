#include "GeneralizedGaussLaguerre1.h"
#include "GaussLegendre.h"

template<int M>
GeneralizedGaussLaguerre1<M>::GeneralizedGaussLaguerre1(std::size_t n, double omega):
    GaussianQuadrature<1, M>(n),
    _omega(omega)
{
    if (n < 1 || n > 10) throw std::domain_error("ERROR: Can only support 1 through 10 integration nodes.");
    if (omega <= 0.0) throw std::invalid_argument("ERROR: Gauss-Laguerre requires that the exponential coefficient be positive");
}

template<int M>
void GeneralizedGaussLaguerre1<M>::setOmega(double omega) const{
    if (omega <= 0.0) throw std::invalid_argument("ERROR: Gauss-Laguerre requires that the exponential coefficient be positive");
    _omega = omega;
}

template<int M>
typename GeneralizedGaussLaguerre1<M>::RangeType GeneralizedGaussLaguerre1<M>::integrate(const Function& f, const IntegrationDomain& D) const{
    RangeType sum = GaussianQuadrature<1, M>::integrate(f, D); // This will give the integral from 0 to infty
    
    // Use a finite-quadrature scheme to compute the integral from 0 to L, and subtract it off
    if (this->eval(D[0]) != 0.0){
        GaussLegendre<1, M> glegendre(n());
        Function f2 = [&f, this](const DomainType& x){ return f(x)*x*std::exp(-_omega*x); };
        IntegrationDomain D2{0.0, D[0]};
        sum -= glegendre.integrate(f2, D2);
    }
    
    return sum;
}

template<int M>
Eigen::VectorXd GeneralizedGaussLaguerre1<M>::getNodes(double, double u) const{
    if (!std::isinf(u) || u < 0) throw std::domain_error("ERROR: Generalized Gauss-Laguerre with a finite upper bound has not been implemented.");
    return _nodes[n()-1] / _omega; // always returns the nodes of the integral from 0 to infty
}

template<int M>
Eigen::VectorXd GeneralizedGaussLaguerre1<M>::getWeights(double, double u) const{
    if (!std::isinf(u) || u < 0) throw std::domain_error("ERROR: Generalized Gauss-Laguerre with a finite upper bound has not been implemented.");
    return _weights[n()-1] / (_omega*_omega); // always returns the weights of the integral from 0 to infty
}