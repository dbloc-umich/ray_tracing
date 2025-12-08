#include "catch2/catch2.hpp"
#include "Ellipsoid.h"
#include "Box.h"
#include "Constants.h"
#include "UnitVector.h"

namespace{
    Ellipsoid sph(Eigen::Vector3d(0.0, 0.0, 0.0), 1, 1, 1); // unit sphere centered at origin
    Ellipsoid el1(Eigen::Vector3d(0.5, 0.0, 0.0), 5, 4, 3); // axis-aligned ellipsoid centered at origin
    Eigen::Matrix3d M{{ 0.08771166, -0.02472796, -0.02203672},
                    {-0.02472796,  0.06735272,  0.00039160},
                    {-0.02203672,  0.00039160,  0.05854673}
                    };

    Ellipsoid el2(Eigen::Vector3d(0.5, 0.0, 0.0), M); // ellipsoid with semi-axes 5, 4, 3
    constexpr double tol = 1e-6;
}

TEST_CASE("ellipsoid_semi_axes_test"){
    REQUIRE(sph.semiAxis(0) == 1);
    REQUIRE(sph.semiAxis(1) == 1);
    REQUIRE(sph.semiAxis(2) == 1);
    REQUIRE(el1.semiAxis(0) == 5);
    REQUIRE(el1.semiAxis(1) == 4);
    REQUIRE(el1.semiAxis(2) == 3);
    REQUIRE(std::abs(1 - el2.semiAxis(0)/5) <= tol);
    REQUIRE(std::abs(1 - el2.semiAxis(1)/4) <= tol);
    REQUIRE(std::abs(1 - el2.semiAxis(2)/3) <= tol);
}

TEST_CASE("ellipsoid_extreme_points_test"){
    REQUIRE(sph.xMin() == -1);
    REQUIRE(sph.xMax() == 1);
    REQUIRE(sph.yMin() == -1);
    REQUIRE(sph.yMax() == 1);
    REQUIRE(sph.zMin() == -1);
    REQUIRE(sph.zMax() == 1);

    REQUIRE(el1.xMin() == -4.5);
    REQUIRE(el1.xMax() == 5.5);
    REQUIRE(el1.yMin() == -4);
    REQUIRE(el1.yMax() == 4);
    REQUIRE(el1.zMin() == -3);
    REQUIRE(el1.zMax() == 3);

    REQUIRE(std::abs(1 + el2.xMin()/3.267660034250652) <= tol);
    REQUIRE(std::abs(1 - el2.xMax()/4.267660034250652) <= tol);
    REQUIRE(std::abs(1 + el2.yMin()/4.091284593459361) <= tol);
    REQUIRE(std::abs(1 - el2.yMax()/4.091284593459361) <= tol);
    REQUIRE(std::abs(1 + el2.zMin()/4.366477784396992) <= tol);
    REQUIRE(std::abs(1 - el2.zMax()/4.366477784396992) <= tol);
}

TEST_CASE("ellipsoid_surface_area_test"){
    double area = 4*mconst::pi;
    REQUIRE(std::abs(1.0 - sph.surfaceArea()/area) <= tol);

    area = 199.455059362;
    REQUIRE(std::abs(1.0 - el1.surfaceArea()/area) <= 1e-4);
    REQUIRE(std::abs(1.0 - el2.surfaceArea()/area) <= 1e-4);
}

TEST_CASE("ellipsoid_volume_test"){
    double volume = 4*mconst::pi/3;
    REQUIRE(std::abs(1.0 - sph.volume()/volume) <= tol);

    volume *= 60;
    REQUIRE(std::abs(1.0 - el1.volume()/volume) <= tol);
    REQUIRE(std::abs(1.0 - el2.volume()/volume) <= tol);
}

TEST_CASE("ellipsoid_point_inside_test"){
    REQUIRE(sph.encloses(Eigen::Vector3d(0.9, 0, 0)));
    REQUIRE(!sph.encloses(Eigen::Vector3d(1.0, 0, 0)));
    REQUIRE(!sph.encloses(Eigen::Vector3d(1.1, 0, 0)));

    REQUIRE(el1.encloses(Eigen::Vector3d(5.4, 0.0, 0.0)));
    REQUIRE(!el1.encloses(Eigen::Vector3d(5.5, 0.0, 0.0)));
    REQUIRE(!el1.encloses(Eigen::Vector3d(5.6, 0.0, 0.0)));

    auto v = el2.principalAxis(0).value();
    REQUIRE(el2.encloses(v*4.9 + Eigen::Vector3d(0.5, 0.0, 0.0)));
    REQUIRE(!el2.encloses(v*5.0 + Eigen::Vector3d(0.5, 0.0, 0.0)));
    REQUIRE(!el2.encloses(v*5.1 + Eigen::Vector3d(0.5, 0.0, 0.0)));
}

TEST_CASE("ellipsoid_point_on_surface_test"){
    REQUIRE(!sph.surfaceContains(Eigen::Vector3d(0.9, 0, 0)));
    REQUIRE(sph.surfaceContains(Eigen::Vector3d(1.0, 0, 0)));
    REQUIRE(!sph.surfaceContains(Eigen::Vector3d(1.1, 0, 0)));

    REQUIRE(!el1.surfaceContains(Eigen::Vector3d(5.4, 0, 0)));
    REQUIRE(el1.surfaceContains(Eigen::Vector3d(5.5, 0, 0)));
    REQUIRE(!el1.surfaceContains(Eigen::Vector3d(5.6, 0, 0)));

    auto v = el2.principalAxis(0).value();
    REQUIRE(!el2.surfaceContains(v*4.9 + Eigen::Vector3d(0.5, 0.0, 0.0)));
    REQUIRE(el2.surfaceContains(v*5.0 + Eigen::Vector3d(0.5, 0.0, 0.0)));
    REQUIRE(!el2.surfaceContains(v*5.1 + Eigen::Vector3d(0.5, 0.0, 0.0)));
}

TEST_CASE("ellipsoid_ellipsoid_intersection_test"){
    REQUIRE(sph.overlaps(el1));
    REQUIRE(sph.overlaps(el2));

    // Barely touching, should return false
    Ellipsoid el3(Eigen::Vector3d(6.0, 0.0, 0.0), 5.0, 4.0, 3.0);
    REQUIRE(!sph.overlaps(el3));

    el3.setSemiAxis(0, 4.0);
    el3.setSemiAxis(1, 5.0);
    REQUIRE(!sph.overlaps(el3));
    REQUIRE(el1.overlaps(el3));
    REQUIRE(el2.overlaps(el3));
}

TEST_CASE("ellipsoid_box_intersection_test"){
    Box box0(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    REQUIRE(el2.overlaps(box0)); // Box center inside ellipsoid

    Box box1(0.0, 0.0, 0.0, 5.0, 5.0, 5.0);
    REQUIRE(el2.overlaps(box1)); // Vertex inside ellipsoid

    Box box2(-5.0, 0, 0, 5.0, 5.0, 5.0);
    REQUIRE(el2.overlaps(box2)); // Edge intersects ellipsoid

    Box box3(-5.0, -5.0, 4.0, 5.0, 5.0, 5.0);
    REQUIRE(el2.overlaps(box3)); // Face intersects ellipsoid

    Box box4(4.0, -5.0, 4.0, 5.0, 5.0, 5.0);
    REQUIRE(!el2.overlaps(box4)); // No intersection
}

TEST_CASE("ellipsoid_distance_to_surface_test"){
    UnitVector3d dir(1.0, 0.0, 0.0);
    // Point inside
    Eigen::Vector3d r(0.0, 0.0, 0.0);
    REQUIRE(std::abs(1 - sph.distanceToSurface(r, dir)/1.0) <= tol);
    REQUIRE(std::abs(1 - el1.distanceToSurface(r, dir)/5.5) <= tol);
    REQUIRE(std::abs(1 - el2.distanceToSurface(r, dir)/3.876535614) <= tol); 
    
    // Point on surface of el2 but outside of el1
    r = Eigen::Vector3d(-1.0, 2.73105203109461, -2.0);
    REQUIRE(std::isnan(sph.distanceToSurface(r, dir)));
    REQUIRE(std::abs(1 - el1.distanceToSurface(r, dir)/0.005090046) <= tol);
    REQUIRE(std::abs(1 - el2.distanceToSurface(r, dir)/3.534932446) <= tol);
    REQUIRE(el2.distanceToSurface(r, -dir) == 0.0);
}

TEST_CASE("ellipsoid_outward_normal_test"){
    REQUIRE(sph.normal(Eigen::Vector3d(1.0, 0.0, 0.0)).isApprox(UnitVector3d(1.0, 0.0, 0.0)));
    REQUIRE(el1.normal(Eigen::Vector3d(5.5, 0.0, 0.0)).isApprox(UnitVector3d(1.0, 0.0, 0.0)));
    //REQUIRE(el2.normal(Eigen::Vector3d(4.2676593302318, 1.3750774494856, 1.40787226197928)).isApprox(UnitVector3d(1.0, 0.0, 0.0)));
}