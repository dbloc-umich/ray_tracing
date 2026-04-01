#include "Shape.h"
#include "Material.h"

#include <iostream>

Shape::Shape(std::shared_ptr<Material> mat):
    _mat(mat)
{}

bool Shape::encloses(const Shape& other) const{
    throw std::runtime_error("ERROR: Containment check cannot be done for non-Box classes.");
    return false;
}

bool Shape::hasProperty(const std::string& name) const noexcept{
    return _mat->hasProperty(name);
}

double Shape::computeProperty(const std::string& name, const std::map<std::string, double>& vars) const{
    return _mat->computeProperty(name, vars);
}

std::ostream& operator<<(std::ostream& os, const Shape& shape){
    return shape.print(os);
}