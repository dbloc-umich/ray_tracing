// #include "catch2/catch2.hpp"
// #include "ParametricSurface.h"
// #include "Constants.h"

// namespace{
//     auto x = [](const Eigen::Vector2d& u){ return std::sqrt(1-u[0]*u[0])*std::cos(u[1]); };
//     auto y = [](const Eigen::Vector2d& u){ return std::sqrt(1-u[0]*u[0])*std::sin(u[1]); };
//     auto z = [](const Eigen::Vector2d& u){ return u[0]; };
//     ParametricSurfaceProperties unitCircle(x, y, z, 0, 1, 0, 2*mconst::pi);

// //     Eigen::Matrix3d M{{ 2.81499928,  0.531864187,  2.44707182},
// //                       { 2.49774303,  2.926071290, -1.39211921},
// //                       { 3.29196841, -2.674925870, -1.03626430}};
// //     Eigen::Vector3d dr{0.5, 0.0, 0.0};
// }
