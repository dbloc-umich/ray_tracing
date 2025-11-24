#include "catch2/catch2.hpp"
#include "Derivative.h"

TEST_CASE("R1-->R1 function"){
    auto f = [](double x){ return x*x; };
    double x = 2.5;
    REQUIRE(std::abs(df(f, x)-5.0) <= 1.0e-8);
    REQUIRE(std::abs(df2(f, x)-2.0) <= 1.0e-8);
}

TEST_CASE("R3-->R1 function"){
    auto f = [](const Eigen::Vector3d& x){ return x[0]*x[0] + x[1]*x[1] + x[2]*x[2]; };
    Eigen::Vector3d x{2, -1, 3};
    Eigen::Vector3d grad_true{4, -2, 6};
    Eigen::Matrix3d hess_true;
    hess_true.setZero();
    hess_true(0,0) = 2.0;
    hess_true(1,1) = 2.0;
    hess_true(2,2) = 2.0;
    REQUIRE((grad(f, x) - grad_true).squaredNorm() <= 1.0e-8);
    REQUIRE((Hessian(f, x)- hess_true).squaredNorm() <= 1.0e-8);
}

TEST_CASE("R3-->R3 function"){
    auto f = [](const Eigen::Vector3d& x){ return Eigen::Vector3d{x[1]*x[1] + x[2]*x[2],
                                                                  -x[0]*x[0] - x[2]*x[2],
                                                                  x[0]*x[0] + x[1]*x[1]};
                                                                  };
    Eigen::Vector3d x{2, -1, 3};
    Eigen::Vector3d curl_true{4, 2, -2};
    Eigen::Matrix3d jac_true;
    jac_true << 0, -2, 6, -4, 0, -6, 4, -2 , 0;
    REQUIRE(std::abs(div(f, x)) <= 1.0e-8);
    REQUIRE((curl(f, x) - curl_true).squaredNorm() <= 1.0e-8);
    REQUIRE((Jacobian(f, x) - jac_true).squaredNorm() <= 1.0e-8);
}