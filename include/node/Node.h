#ifndef NODE_H
#define NODE_H

#include <memory>
#include <vector>

class Shape;
class Node{
    public:
    Node();
    explicit Node(std::unique_ptr<Shape>&& shape); // takes ownership of the pointer

    // Access the shape
    std::unique_ptr<Shape>& shape(){ return _shape; }
    const std::unique_ptr<Shape>& shape() const{ return _shape; }

    // Access the children
    std::unique_ptr<Node>& operator[](std::size_t i){ return _children[i]; }
    const std::unique_ptr<Node>& operator[](std::size_t i) const{ return _children[i]; }

    std::size_t level(){ return _level; }
    void setLevel(std::size_t level){ _level = level; }

    bool contentsOverlap(const Shape& other) const noexcept;
    bool contentsOverlap(const Node& other) const noexcept;

    explicit operator bool() const noexcept{ return static_cast<bool>(_shape); }
    std::size_t size() const noexcept;
    virtual bool empty() const noexcept;
    virtual bool isLeaf() const noexcept{ return _children.empty(); }

    protected:
    std::unique_ptr<Shape> _shape;
    std::vector<std::unique_ptr<Node>> _children;
    std::size_t _level;
};

#endif