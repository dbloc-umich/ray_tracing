#include "FVRobinBC.h"
#include "FVSpatialMesh.h"
#include "NewtonSolver.h"
#include "Eigen/Cholesky"

FVRobinBC::FVRobinBC(const std::function<double(const Eigen::Vector3d&)>& f,
                     const std::function<double(const Eigen::Vector3d&)>& a,
                     const std::function<double(const Eigen::Vector3d&)>& b,
                     std::shared_ptr<FVSpatialMesh> mesh, Eigen::Index surfID,
                     std::shared_ptr<Material> mat, PropVariable var, Prop prop):
    FVBoundaryCondition(mat, var, surfID),
    _f(f),
    _a(a),
    _b(b),
    _mesh(mesh),
    _prop(prop)
{}

double FVRobinBC::computeFlux(double u, const Eigen::Vector3d& r) const{
    // r has to be in the same coordinate system as the problem
    double a = _a(r);
    double b = _b(r) * (_surfID%2 == 0 ? -1 : 1); // lower surface, invert the normal
    double f = _f(r);

    double scale = (_prop == Prop::none) ? 1.0 : _mat->computeProperty(_prop, {{_var, u}});
    double ug; // ghost-cell value
    double dr = (_mesh->ghostLength(_surfID) - _mesh->axis(_surfID/2)[_mesh->axis(_surfID/2).size()-2])/2;
    Eigen::Index i = std::upper_bound(_mesh->axis(0).cbegin(), _mesh->axis(0).cend(), r[0]) - _mesh->axis(0).cbegin() - 1;
    Eigen::Index j = std::upper_bound(_mesh->axis(1).cbegin(), _mesh->axis(1).cend(), r[1]) - _mesh->axis(1).cbegin() - 1;
    Eigen::Index k = std::upper_bound(_mesh->axis(2).cbegin(), _mesh->axis(2).cend(), r[2]) - _mesh->axis(2).cbegin() - 1;

    if (_mesh->isOrthogonal()) ug = (f - u*(a/2 - b/dr))/(a/2 + b/dr);
    else{
        Eigen::Vector3d e0 = _mesh->basisVector(0, r[0], r[1], r[2]);
        Eigen::Vector3d e1 = _mesh->basisVector(1, r[0], r[1], r[2]);
        Eigen::Vector3d e2 = _mesh->basisVector(2, r[0], r[1], r[2]);

        Eigen::Matrix3d J; // Jacobian matrix
        J.col(0) = e0;
        J.col(1) = e1;
        J.col(2) = e2;

        Eigen::Matrix3d G; // Gram matrix
        G(0,0) = e0.dot(e0);
        G(0,1) = e0.dot(e1);
        G(0,2) = e0.dot(e2);
        G(1,0) = G(0,1);
        G(1,1) = e1.dot(e1);
        G(1,2) = e1.dot(e2);
        G(2,0) = G(0,2);
        G(2,1) = G(1,2);
        G(2,2) = e2.dot(e2);
        Eigen::LDLT<Eigen::Ref<Eigen::Matrix3d>> ldlt(G); // LDL^T decomposition

        // These are bounds in the coordinate variables, not necessarily named x, y, or z
        double xL = _mesh->axis(0)[i];
        double xU = _mesh->axis(0)[i+1];
        double yL = _mesh->axis(1)[j];
        double yU = _mesh->axis(1)[j+1];
        double zL = _mesh->axis(2)[k];
        double zU = _mesh->axis(2)[k+1];

        // Function to find the gradient and the boundary value
        auto gradFunc = [&, this](const Eigen::Vector4d& x){
            Eigen::Vector4d rhs;
            Eigen::Vector3d g = x.head(3); // x[0:2] is grad u^b, x[3] is u^g

            Eigen::Vector3d gradRHS;
            gradRHS[0] = (_surfID == 0 || _surfID == 1) ? (x[3]-u)/dr : g.dot(_mesh->cartesian(xU, r[1], r[2]) - _mesh->cartesian(xL, r[1], r[2]))/(xU-xL);
            gradRHS[1] = (_surfID == 2 || _surfID == 3) ? (x[3]-u)/dr : g.dot(_mesh->cartesian(r[0], yU, r[2]) - _mesh->cartesian(r[0], yL, r[2]))/(yU-yL);
            gradRHS[2] = (_surfID == 4 || _surfID == 5) ? (x[3]-u)/dr : g.dot(_mesh->cartesian(r[0], r[1], zU) - _mesh->cartesian(r[0], r[1], zL))/(zU-zL);
            rhs.head(3) = g - J*ldlt.solve(gradRHS);
            
            rhs[3] = a*(x[3]+u)/2 + b*g.dot(_mesh->normalVector(_surfID/2, r[0], r[1], r[2]).value()) - f;
            return rhs;
        };

        NewtonSolver<4> newton(gradFunc); // should be a linear system of equations, but use Newton's regardless
        newton.setFTol(1e-9);
        Eigen::Vector4d grad{u/dr, u/dr, u/dr, u};
        auto status = newton.solve(grad);
        if (status == NLStatus::Success) ug = grad[3];
        throw std::runtime_error("ERROR: Unable to solve for the gradient and ghost-cell value.");
    }

    double scaleg = (_prop == Prop::none) ? 1.0 : _mat->computeProperty(_prop, {{_var, ug}});
    double du = scaleg*ug - scale*u;
    return _mesh->gradientDotN(du/dr, i, j, k, _surfID);
}
