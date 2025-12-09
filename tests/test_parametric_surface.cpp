#include "catch2/catch2.hpp"
#include "ParametricSurface.h"
#include "Constants.h"
#include "Eigen/Eigenvalues"

namespace{
    auto x = [](const Eigen::Vector2d& u){ return std::sin(u[0])*std::cos(u[1]); };
    auto y = [](const Eigen::Vector2d& u){ return std::sin(u[0])*std::sin(u[1]); };
    auto z = [](const Eigen::Vector2d& u){ return std::cos(u[0]); };

    auto xu = [](const Eigen::Vector2d& u){ return std::cos(u[0])*std::cos(u[1]); };
    auto yu = [](const Eigen::Vector2d& u){ return std::cos(u[0])*std::sin(u[1]); };
    auto zu = [](const Eigen::Vector2d& u){ return -std::sin(u[0]); };
    auto xv = [](const Eigen::Vector2d& u){ return -std::sin(u[0])*std::sin(u[1]); };
    auto yv = [](const Eigen::Vector2d& u){ return std::sin(u[0])*std::cos(u[1]); };
    auto zv = [](const Eigen::Vector2d& u){ return 0.0; };

    auto Hx = [](const Eigen::Vector2d& u){
        Eigen::Matrix2d H;
        H(0,0) = -std::sin(u[0])*std::cos(u[1]);
        H(0,1) = -std::cos(u[0])*std::sin(u[1]);
        H(1,0) = H(0,1);
        H(1,1) = H(0,0);
        return H;
    };

    auto Hy = [](const Eigen::Vector2d& u){
        Eigen::Matrix2d H;
        H(0,0) = -std::sin(u[0])*std::sin(u[1]);
        H(0,1) = -std::cos(u[0])*std::cos(u[1]);
        H(1,0) = H(0,1);
        H(1,1) = H(0,0);
        return H;
    };

    auto Hz = [](const Eigen::Vector2d& u){
        Eigen::Matrix2d H = Eigen::Matrix2d::Zero();
        H(0,0) = -std::cos(u[0]);
        return H;
    };

    auto base = std::make_shared<ParametricSurfaceProperties>(x, y, z, 0, mconst::pi, 0, 2*mconst::pi, 2, 4, xu, yu, zu, xv, yv, zv, Hx, Hy, Hz);

    Eigen::Matrix3d M{{ 0.08771166, -0.02472796, -0.02203672},
                      {-0.02472796,  0.06735272,  0.00039160},
                      {-0.02203672,  0.00039160,  0.05854673}}; // same matrix used in test_ellipsoid;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(M);
    Eigen::Matrix3d Q = es.eigenvectors();
    Eigen::Array3d lambda = es.eigenvalues();
    Eigen::Matrix3d A = Q * Eigen::DiagonalMatrix<double, 3>(1.0 / Eigen::sqrt(lambda)) * Q.transpose();
    Eigen::Vector3d dr(0.5, 0.0, 0.0);
    ParametricSurface ellipsoid(base, A, dr);
    constexpr double tol = 1e-6;
}

TEST_CASE("parametric_extreme_points_test"){
    REQUIRE(std::abs(1 + ellipsoid.xMin()/3.267660034250652) <= tol);
    REQUIRE(std::abs(1 - ellipsoid.xMax()/4.267660034250652) <= tol);
    REQUIRE(std::abs(1 + ellipsoid.yMin()/4.091284593459361) <= tol);
    REQUIRE(std::abs(1 - ellipsoid.yMax()/4.091284593459361) <= tol);
    REQUIRE(std::abs(1 + ellipsoid.zMin()/4.366477784396992) <= tol);
    REQUIRE(std::abs(1 - ellipsoid.zMax()/4.366477784396992) <= tol);
}
