#include "InterpolatedProperty.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

InterpolatedProperty::InterpolatedProperty(std::string file)
{
    loadFile(file);
    for (auto it = _grid.begin(); it != std::prev(_grid.end()); ++it){
        if (*it >= *std::next(it))
            throw std::domain_error("ERROR: Independent variable values are not sorted in strictly ascending order.");
    }
    if (_grid.size() != _val.size())
        throw std::invalid_argument("ERROR: Dimension mismatch between independent and dependent variables");
    if (!std::is_sorted(_grid.begin(), _grid.end())){
        // sort vectors here
    }
}

double InterpolatedProperty::compute(const std::vector<double>& vars) const{
    // Can only deal with 1D interpolation for now
    double x = vars[0];
    if (x < _grid[0] || x > _grid.back())
        std::cerr << "WARNING: Variable value is outside of the interpolating range." << std::endl;
    if (x <= _grid[0]) return _val[0];
    if (x >= _grid.back()) return _val.back();
    std::size_t ind = 0;
    for (auto it = _grid.crbegin(); it != _grid.crend(); ++it){
        if (*it <= x){
            ind = std::distance(_grid.begin(), it.base()) - 1;
            break;
        }
    }
    return _val[ind] + (_val[ind+1]-_val[ind]) / (_grid[ind+1]-_grid[ind]) * (x - _grid[ind]); 
}

void InterpolatedProperty::loadFile(std::string file){
    std::ifstream input(file);
    if (!input.is_open()){
        auto msg = "ERROR: Could not open file \"" + file + "\"."; 
        throw std::runtime_error(msg);
    }

    std::string line;
    if (std::getline(input, line)){
        std::istringstream iss(line);
        double num;
        while (iss >> num) _grid.push_back(num);
    }

    _val.reserve(_grid.size());
    if (std::getline(input, line)){
        std::istringstream iss(line);
        double num;
        while (iss >> num) _val.push_back(num);
    }
}