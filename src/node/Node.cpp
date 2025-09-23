#include "Node.h"
#include "Shape.h"

Node::Node():
    _shape(),
    _children{},
    _level(0)
{}

explicit Node::Node(std::unique_ptr<Shape>&& shape):
    _shape(std::move(shape)),
    _children{},
    _level(0)
{}

std::size_t Node::size() const noexcept{
    std::size_t sz = 0;
    for (const auto& node: _children){
        if (node) sz++;
    }
    return sz;
}

bool Node::empty() const noexcept{
    for (const auto& node: _children){
        if (node) return false;
    }
    return true;
}