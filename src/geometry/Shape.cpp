#include "Shape.h"
#include "Material.h"

#include <iostream>

Shape::Shape(std::shared_ptr<Material> mat):
    _mat(mat)
{}

bool Shape::hasProperty(const Prop& name) const noexcept{
    return _mat->hasProperty(name);
}

double Shape::computeProperty(const Prop& name, const std::vector<double>& vars) const{
    return _mat->computeProperty(name, vars);
}

std::ostream& operator<<(std::ostream& os, const Shape& shape){
    return shape.print(os);
}