#ifndef NODE_H
#define NODE_H

#include <memory>
#include <vector>

class Shape;
class Node{
    public:
    explicit Node(std::size_t N = 2);
    explicit Node(std::unique_ptr<Shape> shape, std::size_t N = 2); // takes ownership of the pointer

    // Access the children
    std::unique_ptr<Node>& operator[](std::size_t i){ return _children[i]; }
    const std::unique_ptr<Node>& operator[](std::size_t i) const{ return _children[i]; }

    std::size_t level(){ return _level; }
    void setLevel(std::size_t level){ _level = level; }

    virtual bool leavesOverlap(const Shape& other) const noexcept;
    virtual bool leavesOverlap(const Node& other) const noexcept;

    std::size_t size() const noexcept;
    std::size_t max_size() const noexcept{ return _children.size(); }
    virtual bool empty() const noexcept;
    virtual bool isLeaf() const noexcept{ return empty(); }

    // Modifier and observer functions, from std::unique_ptr
    Shape* release() noexcept{ return _shape.release(); }
    void reset(Shape* ptr = nullptr) noexcept{ _shape.reset(ptr); }
    void swap(std::unique_ptr<Shape>& other) noexcept{ _shape.swap(other); }
    Shape* get() const noexcept{ return _shape.get(); }
    explicit operator bool() const noexcept{ return static_cast<bool>(_shape); }
    Shape& operator*() const noexcept{ return *_shape; }
    Shape* operator->() const noexcept{ return _shape.get(); }

    friend std::ostream& operator<<(std::ostream& os, const Node& node){ return node.print(os); }

    protected:
    std::unique_ptr<Shape> _shape;
    std::vector<std::unique_ptr<Node>> _children;
    std::size_t _level;

    std::ostream& printTabs(std::ostream& os, std::size_t n) const noexcept;
    virtual std::ostream& print(std::ostream& os) const noexcept;
};

#endif