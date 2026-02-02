#include "GaussKronrod.h"

template<int N, int M>
typename GaussKronrod<N, M>::RangeType GaussKronrod<N, M>::integrate(const Function& f, const IntegrationDomain& D) const{
    if (!this->validIntegrationDomain(D)) throw std::invalid_argument("ERROR: Invalid integration domain.");
    return adaptiveIntegral(f, D, _tol, 0);
}

template<int N, int M>
Eigen::VectorXd GaussKronrod<N, M>::getNodes(double l, double u) const{
    Eigen::ArrayXd x(n());
    x.head(this->_n) = _nodes[this->_n-1];
    x.tail(this->_n+1) = _Knodes[this->_n-1];
    return x * (u-l)/2 + (u+l)/2;    
}

template<int N, int M>
Eigen::VectorXd GaussKronrod<N, M>::getWeights(double l, double u) const{
    Eigen::ArrayXd w(n());
    w.head(this->_n) = _Lweights[this->_n-1];
    w.tail(this->_n+1) = _Kweights[this->_n-1];
    return w * (u-l)/2;   
}

template<int N, int M>
double GaussKronrod<N, M>::error(const RangeType& GL, const RangeType& GK) const noexcept{
    // Empirical error estimation using Gauss-Legendre (GL)
    if constexpr (M == 1) return (GK == 0.0) ? std::abs(GL) : std::abs((GL-GK)/GK);
    else{
        RangeType diff = GL-GK;
        for (Eigen::Index i = 0; i < diff.size(); i++) diff[i] /= (GK[i] == 0) ? 1 : GK[i];
        return diff.norm();
    }
}

template<int N, int M>
typename GaussKronrod<N, M>::IntegrationBound GaussKronrod<N, M>::midpoint(const IntegrationBound& l, const IntegrationBound& u) const{
    // Midpoint integration bound, used in adaptive steps
    if (std::get_if<double>(&l)){
        if (std::get_if<double>(&u)) return (this->eval(l) + this->eval(u)) / 2;
        if (std::get_if<IFunc1D>(&u)) return [&l, &u, this](double x) { return (this->eval(l) + this->eval(u,x)) / 2; };
        return [&l, &u, this](const Eigen::VectorXd& x) { return (this->eval(l) + this->eval(u,x)) / 2; };
    } else if (std::get_if<IFunc1D>(&l)) return [&l, &u, this](double x) { return (this->eval(l,x) + this->eval(u,x)) / 2; };
    return [&l, &u, this](const Eigen::VectorXd& x) { return (this->eval(l,x) + this->eval(u,x)) / 2; };
}

template<int N, int M>
void GaussKronrod<N, M>::recursiveSplitDomain(std::vector<IntegrationDomain>& newDomains, const IntegrationDomain& oldDomain,
                              const std::vector<IntegrationBound>& midpoints, std::size_t d, std::size_t start, std::size_t stop) const noexcept{
    if (d == midpoints.size()) return;
    auto pConst = std::get_if<double>(&midpoints[d]);
    if (pConst && std::isnan(*pConst)){
        // No split, midpoint stored as a nan
        for (std::size_t i = start; i < stop; i++){
            newDomains[i][2*d] = oldDomain[2*d];
            newDomains[i][2*d+1] = oldDomain[2*d+1];
        }
        recursiveSplitDomain(newDomains, oldDomain, midpoints, d+1, start, stop);
    } else{
        std::size_t diff = (stop-start)/2;
        for (std::size_t i = start; i < stop; i++){
            newDomains[i][2*d] = (i < start+diff) ? oldDomain[2*d] : midpoints[d];
            newDomains[i][2*d+1] = (i < start+diff) ? midpoints[d] : oldDomain[2*d+1];
        }
        recursiveSplitDomain(newDomains, oldDomain, midpoints, d+1, start, start+diff);
        recursiveSplitDomain(newDomains, oldDomain, midpoints, d+1, start+diff, stop);
    }
}

template<int N, int M>
typename GaussKronrod<N, M>::RangeType GaussKronrod<N, M>::adaptiveIntegral(const Function& f, const IntegrationDomain& D, double tol, std::size_t d) const{
    RangeType sum = GaussianQuadrature<N,M>::integrate(f, D);
    GaussLegendre<N,M> GL(this->_n);
    // would be nice to reuse repeated nodes, but in higher dimensions it doesn't help much
    double err = error(GL.integrate(f, D), sum);

    // Not yet converged, do adaptive quadrature
    if (err > tol){
        if constexpr (N == 1){
            double L = std::get<double>(D[0]);
            double U = std::get<double>(D[1]);
            std::size_t m = std::ceil(std::pow(err/tol, 1.0/(3*this->_n+1))); // number of subintervals
            double dx = (U-L)/m;
            sum *= 0; // reset the accumulator
            for (std::size_t i = 0; i < m; i++) sum += adaptiveIntegral(f, {L+i*dx, L+(i+1)*dx}, tol/m, d);
        } else{
            constexpr std::size_t dim = (N != Eigen::Dynamic) ? N : D.size();
            std::size_t numSplits = std::ceil(std::log2(err/tol) / (3*this->_n+1)); // number of dimensions to split
            if (numSplits > dim) numSplits = dim;
            std::vector<IntegrationBound> midpoints(dim, std::numeric_limits<double>::quiet_NaN()); // stores the midpoints of the split dimensions
            for (std::size_t ns = 0; ns < numSplits; ns++){
                midpoints[d] = midpoint(D[2*d], D[2*d+1]);
                d = (d+1) % dim;
            }

            std::size_t m = std::pow(2, numSplits); // number of subdomains
            std::vector<IntegrationDomain> newDomains(m, IntegrationDomain(2*dim));
            recursiveSplitDomain(newDomains, D, midpoints, 0, 0, m);
            sum *= 0.0; // reset the accumulator
            for (auto domain: newDomains) sum += adaptiveIntegral(f, domain, tol/m, d);
        }
    }
        
    return sum;
}
