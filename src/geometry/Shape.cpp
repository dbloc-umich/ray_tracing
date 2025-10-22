#include "Shape.h"
#include "Material.h"

Shape::Shape(const std::shared_ptr<Material>& mat):
    _mat(mat)
{}

double Shape::getProp(std::string name, const std::vector<double>& vars) const{
    return _mat->computeProperty(name, vars);
}

std::ostream& operator<<(std::ostream& os, const Shape& shape){
    return shape.print(os);
}