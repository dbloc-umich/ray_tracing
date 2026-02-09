// #include "InterpolatedProperty.h"

// #include <algorithm>
// #include <stdexcept>

// InterpolatedProperty::InterpolatedProperty(const Eigen::ArrayXd& grid, const Eigen::ArrayXd& val):
//     _grid(grid),
//     _val(val)
// {
//     if (std::is_sorted(_grid.cbegin(), _grid.cend(), [](double a, double b){ return a < b; }))
//         throw std::domain_error("ERROR: Independent variable values are not sorted in strictly ascending order.");
//     if (_grid.size() != _val.size())
//         throw std::invalid_argument("ERROR: Dimension mismatch between independent and dependent variables");
// }

// double InterpolatedProperty::compute(const Eigen::ArrayXd& vars) const{
//     // Can only deal with 1D linear interpolation for now
//     double x = vars[0];
//     // if (x < _grid[0] || x > _grid[_grid.size()-1])
//     //     std::cerr << "WARNING: Variable value is outside of the interpolating range." << std::endl;
//     if (x <= _grid[0]) return _val[0];
//     if (x >= _grid[_grid.size()-1]) return _val[_grid.size()-1];
    
//     Eigen::Index i = std::upper_bound(_grid.cbegin(), _grid.cend(), x) - _grid.cbegin();
//     return _val[i-1] + (_val[i]-_val[i-1]) / (_grid[i]-_grid[i-1]) * (x - _grid[i-1]); 
// }