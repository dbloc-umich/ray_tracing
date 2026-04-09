/**
 * Using generalized Gauss-Laguerre Quadrature to evaluate the integral of x*f(x)exp(-omega*x) over a semi-infinite domain [l, infty].
 * Only alpha = 1 and N = 1 (single integral) cases have been implemented
 * Attempting to use Gauss-Laguerre to integrate a general function g(x) by setting f(x) = x*g(x)exp(omega*x) may result in instabilities.
**/

#ifndef GENERALIZED_GAUSS_LAGUERRE_H
#define GENERALIZED_GAUSS_LAGUERRE_H

#include "GaussianQuadrature.h"

class GeneralizedGaussLaguerre1Base{
    // Dummy class to store the nodes and weights
    public:
    virtual ~GeneralizedGaussLaguerre1Base() = default;
    inline static const std::vector<Eigen::ArrayXd> _nodes =
        {{{2.0000000000000000}}, // starting at n = 1
         {{1.2679491924311228, 4.732050807568877}},
         {{0.9358222275240875, 3.3054072893322790, 7.758770483143635}},
         {{0.7432919279814324, 2.5716350076462793, 5.7311787516890975, 10.953894312683188}},
         {{0.6170308532782699, 2.1129659585785228, 4.6108331510175310,  8.399066971204840, 14.260103065920832}}, // n = 5
         {{0.5276681217111288, 1.7962998096434084, 3.8766415204769125,  6.918816566704721, 11.234610429083117, 17.645963552380710}},
         {{0.4610242198049950, 1.5635861896542630, 3.3520505025367333,  5.916297249020539,  9.420699383021590, 14.194165548007465, 21.092176907954414}},
         {{0.4093835732031853, 1.3849631848031398, 2.9562545561688620,  5.181943101040074,  8.161709688145818, 12.070055126837152, 17.249735526149000, 24.585955243652784}},
         {{0.3681784529417420, 1.2433579621404816, 2.6460338413842077,  4.616882514635060,  7.221786539396572, 10.567320807741863, 14.835914515260932, 20.382181985449247, 28.118343381049886}},
         {{0.3345286763247511, 1.1282533558766334, 2.3958699247473043,  4.166840987928769,  6.487353031380815,  9.428354813335615, 13.101723580367802, 17.696487566846226, 23.577787088360154, 31.682800974831935}}
        }; // n = 10 is the highest order
    
    inline static const std::vector<Eigen::ArrayXd> _weights =
        {{{1.0000000000000000}}, // starting at n = 1
         {{0.7886751345948129, 0.2113248654051872}},
         {{0.5886814810396593, 0.3912160592223104, 0.0201024597380305}},
         {{0.4468705932187766, 0.4776357723638678, 0.0741777847310521, 0.0013158496863032}},
         {{0.3480145400233489, 0.5022806741324928, 0.1409159194944728, 0.0087198930260100, 6.8973323585640127e-05}}, // n = 5
         {{0.2776501420298750, 0.4939105830503542, 0.2030042967437300, 0.0246688203691898, 7.6304276746352873e-04, 3.1150393875275513e-06}},
         {{0.2262105963831891, 0.4697087077412667, 0.2531225162324003, 0.0478031088842370, 3.1004597560253984e-03, 5.4484668893413588e-05, 1.2633398815362801e-07}},
         {{0.1876325414057236, 0.4389853607311424, 0.2899960707813138, 0.0751413846166972, 7.9326466487073376e-03, 3.0864213681330502e-04, 3.3489582097970791e-06, 4.7213928231932295e-09}},
         {{0.1580309235007815, 0.4065825794555731, 0.3149824836650979, 0.1037755199516433, 1.5577637724278439e-02, 1.0248719698507575e-03, 2.5800211387172334e-05, 1.8335595501724325e-07, 1.6543294559699223e-10}},
         {{0.1348565545313764, 0.3749243411054594, 0.3302471308042169, 0.1315288846413067, 2.5842064211717062e-02, 2.4896457958959364e-03, 1.0948847642669560e-04, 1.8812754141816992e-06, 9.1526850285756175e-09, 5.5013034961787156e-12}}
        }; // n = 10 is the highest order
};

template<int M = 1>
class GeneralizedGaussLaguerre1 : public GaussianQuadrature<1, M>, GeneralizedGaussLaguerre1Base{
    public:
    using typename Quadrature<1, M>::DomainType;
    using typename Quadrature<1, M>::RangeType;
    using typename Quadrature<1, M>::Function;
    using typename Quadrature<1, M>::IntegrationDomain;

    GeneralizedGaussLaguerre1(std::size_t n, double omega=1.0);
    double omega() const noexcept{ return _omega; }
    void setOmega(double) const;
    RangeType integrate(const Function& f, const IntegrationDomain& D) const override;

    protected:
    mutable double _omega;
    std::size_t n() const noexcept override{ return this->_n; }
    Eigen::VectorXd getNodes(double l, double u) const override;
    Eigen::VectorXd getWeights(double l, double u) const override;
};

#include "GeneralizedGaussLaguerre1.tpp"

#endif