#include "FVSpatialMesh.h"
#include <algorithm>

FVSpatialMesh::FVSpatialMesh(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z):
    _axes{x, y, z}
{
    auto comp = [](double a, double b){ return a < b; };
    for (std::size_t i = 0; i < 3; i++){
        if (!std::is_sorted(_axes[i].cbegin(), _axes[i].cend(), comp))
            throw std::invalid_argument("ERROR: Axis boundary values are not sorted in strictly ascending order.");
        if (_axes[i].size() < 2)
            throw std::invalid_argument("ERROR: Axis has too few entries.");
    }
}

FVSpatialMesh::FVSpatialMesh(const std::array<Eigen::ArrayXd, 3>& axes):
    _axes(axes)
{
    auto comp = [](double a, double b){ return a < b; };
    for (std::size_t i = 0; i < 3; i++){
        if (!std::is_sorted(_axes[i].cbegin(), _axes[i].cend(), comp))
            throw std::invalid_argument("ERROR: Axis boundary values are not sorted in strictly ascending order.");
        if (_axes[i].size() < 2)
            throw std::invalid_argument("ERROR: Axis has too few entries.");
    }
}

