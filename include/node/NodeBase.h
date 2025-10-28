/**
 * This is a base template for Node classes, which contains a pointer to a Shape object and multiple pointers to its child Nodes.
 * The child nodes must be of the same class as the parent Node (i.e., the Tree design is homogeneous), leading to the CRTP Node design
 * The simplest Node class does not deviate from this behavior. However, more complex Node designs may include:
 *  An octree Node that contains its own leaves in addition to the 8 octants
 *  An octree Node that shares some of its leaves with other Nodes
 *  A bidirectional Node that has a pointer back to its parent
**/

#ifndef NODE_BASE_H
#define NODE_BASE_H

#include <algorithm>
#include <array>
#include <iostream>
#include <memory>

class Shape;

template<class T, std::size_t N = 2>
class NodeBase{
    public:
    explicit NodeBase(Shape* shape = nullptr): _shape(shape), _children(), _level(0) {}
    explicit NodeBase(std::unique_ptr<Shape> shape = nullptr): _shape(std::move(shape)), _children(), _level(0) {}
    NodeBase(const NodeBase<T,N>&) = delete;
    NodeBase(NodeBase<T,N>&&) = default;
    virtual ~NodeBase() = default;
    NodeBase<T,N>& operator=(const NodeBase<T,N>&) = delete;
    NodeBase<T,N>& operator=(NodeBase<T,N>&&) = default;

    // Access the children
    std::unique_ptr<T>& operator[](std::size_t i){ return _children[i]; }
    const std::unique_ptr<T>& operator[](std::size_t i) const{ return _children[i]; }

    std::size_t level() const noexcept{ return _level; }
    void setLevel(std::size_t level) noexcept{ _level = level; }

    virtual bool leavesOverlap(const Shape& other) const noexcept = 0;
    virtual bool leavesOverlap(const T& other) const noexcept = 0;

    std::size_t size() const noexcept{
        return std::count_if(_children.cbegin(), _children.cend(), [](auto& node){ return bool(node); }); }
    constexpr std::size_t max_size() const noexcept{ return N; }
    virtual bool empty() const noexcept{
        return std::find_if(_children.cbegin(), _children.cend(), [](auto& node){ return bool(node); }) == _children.cend(); }
    virtual bool isLeaf() const noexcept{ return empty(); }

    // Modifier and observer functions, from std::unique_ptr
    Shape* release() noexcept{ return _shape.release(); }
    void reset(Shape* ptr = nullptr) noexcept{ _shape.reset(ptr); }
    void swap(std::unique_ptr<Shape>& other) noexcept{ _shape.swap(other); }
    Shape* get() const noexcept{ return _shape.get(); }
    explicit operator bool() const noexcept{ return static_cast<bool>(_shape); }
    Shape& operator*() const noexcept{ return *_shape; }
    Shape* operator->() const noexcept{ return _shape.get(); }

    friend std::ostream& operator<<(std::ostream& os, const NodeBase<T, N>& node){ return node.print(os); }

    protected:
    std::unique_ptr<Shape> _shape;
    std::array<std::unique_ptr<T>, N> _children;
    std::size_t _level;

    std::ostream& printTabs(std::ostream& os, std::size_t n) const noexcept;
    virtual std::ostream& print(std::ostream& os) const noexcept = 0;
};

template<class T, std::size_t N>
std::ostream& NodeBase<T, N>::printTabs(std::ostream& os, std::size_t n) const noexcept{
    for (unsigned i = 0; i < n; i++) os << "  ";
    return os;
}

#endif