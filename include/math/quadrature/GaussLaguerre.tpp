#include "GaussLaguerre.h"

template<int N, int M>
GaussLaguerre<N, M>::GaussLaguerre(std::size_t n, double omega):
    GaussianQuadrature<N, M>(n),
    _omega(omega)
{
    if (n < 1 || n > 10) throw std::domain_error("ERROR: Can only support 1 through 10 integration nodes.");
    if (omega <= 0.0) throw std::invalid_argument("ERROR: Gauss-Laguerre requires that the exponential coefficient be positive");
}

template<int N, int M>
void GaussLaguerre<N, M>::setOmega(double omega) const{
    if (omega <= 0.0) throw std::invalid_argument("ERROR: Gauss-Laguerre requires that the exponential coefficient be positive");
    _omega = omega;
}

template<int N, int M>
Eigen::VectorXd GaussLaguerre<N, M>::getNodes(double l, double u) const{
    Eigen::ArrayXd x;
    if (std::isinf(u) && u > 0){
        x.resize(n());
        x = _nodes[n()-1] / _omega + l;
    } else{
        x.resize(2*n());
        x.head(n()) = _nodes[n()-1] / _omega + l;
        x.tail(n()) = _nodes[n()-1] / _omega + u;
    }
    return x;
}

template<int N, int M>
Eigen::VectorXd GaussLaguerre<N, M>::getWeights(double l, double u) const{
    Eigen::ArrayXd w;
    if (std::isinf(u) && u > 0){
        w.resize(n());
        w = _weights[n()-1] * std::exp(-_omega*l)/_omega;
    } else{
        w.resize(2*n());
        w.head(n()) = _weights[n()-1] * std::exp(-_omega*l)/_omega;
        w.tail(n()) = -_weights[n()-1] * std::exp(-_omega*u)/_omega;
    }
    return w;
}
