// #include "ParametricSurfaceProperties.h"
// #include "NonlinearSolver.h"
// #include <algorithm>

// ParametricSurfaceProperties::ParametricSurfaceProperties(const function_type& x, const function_type& y, const function_type& z,
//                                                         double uMin, double uMax, double vMin, double vMax,
//                                                         //unsigned short uSymmetry, unsigned short vSymmetry,
//                                                         const function_type& xu, const function_type& yu, const function_type& zu,
//                                                         const function_type& xv, const function_type& yv, const function_type& zv,
//                                                         const Hfunction_type& Hx, const Hfunction_type& Hy, const Hfunction_type& Hz):
//     FunctionShapeProperties(),
//     _x(x), _y(y), _z(z),
//     _u0(std::min(uMin, uMax)), _u1(std::max(uMin, uMax)),
//     _v0(std::min(vMin, vMax)), _v1(std::max(vMin, vMax)),
//     //_uSym(uSymmetry), _vSym(vSymmetry),
//     _xu(xu ? xu : [this](auto& u){ return df(_x, u, 0); }),
//     _yu(yu ? yu : [this](auto& u){ return df(_y, u, 0); }),
//     _zu(zu ? zu : [this](auto& u){ return df(_z, u, 0); }),
//     _xv(xv ? xv : [this](auto& u){ return df(_x, u, 1); }),
//     _yv(yv ? yv : [this](auto& u){ return df(_y, u, 1); }),
//     _zv(zv ? zv : [this](auto& u){ return df(_z, u, 1); }),
//     _Hx(Hx ? Hx : [this](auto& u){ return Hessian(_x, u); }),
//     _Hy(Hy ? Hy : [this](auto& u){ return Hessian(_y, u); }),
//     _Hz(Hz ? Hz : [this](auto& u){ return Hessian(_z, u); })
// {
//     if (!_x || !_y || !_z)
//         throw std::invalid_argument("ERROR: Invalid parametrization of at least one spatial variable.");
//     if (_u0 == _u1 || _v0 == _v1)
//         throw std::invalid_argument("ERROR: Repeated bound in at least one parametric variable.");
//     // if (_uSym == 0 || _vSym == 0)
//     //     throw std::invalid_argument("ERROR: Invalid _uSym or _vSym values; must be a positive integer.");
//     checkClosedSurface();
//     computeExtrema();
//     computeSurfaceArea();
//     if (_isClosedSurface) computeVolume();
//     else _volume = std::numeric_limits<double>::quiet_NaN();
// }

// void ParametricSurfaceProperties::checkClosedSurface(){
//     // Sample points
//     constexpr Eigen::Index N = 11;
//     std::array<Eigen::VectorXd, 2> samplePoints;
//     samplePoints[0].resize(N);
//     samplePoints[1].resize(N);
//     for (Eigen::Index i = 0; i < N; i++){
//         samplePoints[0][i] = _v0 + (_v1-_v0)/(N-1)*i;
//         samplePoints[1][i] = _u0 + (_u1-_u0)/(N-1)*i;
//     }

//     for (Eigen::Index i = 0; i < 3; i++){
//         for (Eigen::Index j = 0; j < 2; j++){
//             if (!isPeriodic(i, j, samplePoints[j]) && !isDegenerate(i, j, samplePoints[j])){
//                 _isClosedSurface = false;
//                 return;
//             }
//         }
//     }
//     _isClosedSurface = true;
// };

// void ParametricSurfaceProperties::computeExtrema(){
//     for (std::size_t i = 0; i < 3; i++){
//         // Check the corners first      
//         auto x = (i == 0 ? _x : ((i == 1) ? _y : _z)); // Spatial variable to optimize   
//         double xMin = std::numeric_limits<double>::max();
//         double xMax = std::numeric_limits<double>::lowest();
//         for (double u: {_u0, _u1}){
//             for (double v: {_v0, _v1}){
//                 double val = x({u, v});
//                 if (val < xMin) xMin = val;
//                 if (val > xMax) xMax = val;
//             }
//         }

//         auto grad = [this, &i](const auto& u){
//             if (i == 0) return gradx(u);
//             if (i == 1) return grady(u);
//             return gradz(u);
//         };
//         auto& H = (i == 0) ? _Hx : ((i == 1) ? _Hy : _Hz);

//         // Sample values of [u, v] to find all critical points
//         std::vector<Eigen::Vector2d> critPts;
//         constexpr std::size_t Nu = 6;
//         constexpr std::size_t Nv = 6;

//         // Using Newton's method to find all possible critical points
//         for (std::size_t i = 0; i < Nu; i++){
//             for (std::size_t j = 0; j < Nv; j++){
//                 Eigen::Vector2d crit{_u0+(_u1-_u0)/(Nu-1)*i, _v0+(_v1-_v0)/(Nv-1)*j};
//                 try{
//                     newton(grad, H, crit);
//                     if (crit[0] >= _u0 && crit[0] <= _u1 && crit[1] >= _v0 && crit[1] <= _v1 ){
//                         double val = x(crit);
//                         if (val < xMin) xMin = val;
//                         if (val > xMax) xMax = val;
//                     }
//                 } catch(const std::runtime_error& ex){};
//             }
//         }

//         if (i == 0){
//             _xMin = xMin;
//             _xMax = xMax;
//         } else if (i == 1){
//             _yMin = xMin;
//             _yMax = xMax;
//         } else{
//             _zMin = xMin;
//             _zMax = xMax;
//         }
//     }
// }

// void ParametricSurfaceProperties::computeSurfaceArea(){
//     auto dS = [this](const Eigen::Vector2d& u){ return ru(u).cross(rv(u)).norm(); };
    
//     // Gauss-Legendre integration with 3 points
//     _surfaceArea = 0.0;
//     double uub = _u1; //_u0 + (_u1-_u0)/_uSym; // upper bound in u
//     double vub = _v1; //_v0 + (_v1-_v0)/_vSym; // upper bound in v
//     const Eigen::Array3d roots{-std::sqrt(3.0/5), 0.0, std::sqrt(3.0/5)};
//     const Eigen::Array3d weights{5.0/9, 8.0/9, 5.0/9};
    
//     for (Eigen::Index i = 0; i < 3; i++){
//         double u = ((uub-_u0)*roots[i] + uub+_u0)/2;
//         double wu = (uub-_u0)/2 * weights[i];
//         for (Eigen::Index j = 0; j < 3; j++){
//             double v = ((vub-_v0)*roots[j] + vub+_v0)/2;
//             double wv = (vub-_v0)/2 * weights[j];
//             _surfaceArea += wu * wv * dS({u,v});
//         }
//     }
//     //_surfaceArea *= _uSym*_vSym;
// }

// void ParametricSurfaceProperties::computeVolume(){
//     auto dV = [this](const Eigen::Vector2d& u){ return _z(u)*(_xu(u)*_yv(u) - _yu(u)*_xv(u)); };

//     // Gauss-Legendre integration with 3 points
//     _volume = 0.0;
//     const Eigen::Array3d roots{-std::sqrt(3.0/5), 0.0, std::sqrt(3.0/5)};
//     const Eigen::Array3d weights{5.0/9, 8.0/9, 5.0/9};
    
//     for (Eigen::Index i = 0; i < 3; i++){
//         double u = ((_u1-_u0)*roots[i] + _u1+_u0)/2;
//         double wu = (_u1-_u0)/2 * weights[i];
//         for (Eigen::Index j = 0; j < 3; j++){
//             double v = ((_v1-_v0)*roots[j] + _v1+_v0)/2;
//             double wv = (_v1-_v0)/2 * weights[j];
//             _volume += wu * wv * dV({u,v});
//         }
//     }
//     _volume = std::abs(_volume);
// };

// bool ParametricSurfaceProperties::isPeriodic(Eigen::Index i, Eigen::Index j, const Eigen::VectorXd& arr){
//     // Check if a coordinate variable is periodic (e.g. x(u0,v) == x(u1,v) for all v)
//     auto& x = (i == 0 ? _x : ((i == 1) ? _y : _z)); // Spatial variable to operate on
//     auto f = [this, j, &x](double t){ // Helper function to check if the function is periodic at each sample point
//         double x0 = j == 0 ? x({_u0, t}) : x({t, _v0});
//         double x1 = j == 0 ? x({_u1, t}) : x({t, _v1});
//         return std::abs(x0-x1) <= std::max(1.0, std::max(x0,x1))*FunctionShapeProperties::eps;
//     };

//     for (auto t: arr){
//         if (!f(t)) return false;
//     }
//     return true;
// }

// bool ParametricSurfaceProperties::isDegenerate(Eigen::Index i, Eigen::Index j, const Eigen::VectorXd& arr){
//     // Check if a coordinate variable is degenerate to a point (e.g. x(u0,v) and x(u1,v) are not functions v)
//     auto& dx = (j == 0 ? (i == 0 ? _xv : (i == 1) ? _yv : _zv)
//                 : (i == 0 ? _xu : (i == 1) ? _yu : _zu)); // Spatial variable to operate on
//     // Helper functions to check if the function is degenerate at each endpoint
//     auto f = [this, j, &dx](double t){
//         double dx0 = (j == 0) ? dx({_u0, t}) : dx({t, _v0});
//         return std::abs(dx0) <= FunctionShapeProperties::eps;
//     };
//     auto g = [this, j, &dx](double t){
//         double dx0 = (j == 0) ? dx({_u1, t}) : dx({t, _v1});
//         return std::abs(dx0) <= FunctionShapeProperties::eps;
//     };
//     for (auto t: arr){
//         if (!f(t) || !g(t)) return false;
//     }
//     return true;
// }