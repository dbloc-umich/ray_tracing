#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H
#include "Box.h"
#include<array>
#include <memory>
#include <type_traits>
#include <vector>

using Node = std::unique_ptr<Shape>; // may be redefined as a class later to include bidirectional nodes

class BoundingBox: public Box{
    public:
    /* To be worked on later
    template <typename... Args,
              typename = std::enable_if_t<!std::is_base_of<BoundingBox, std::decay_t<Args>...>::value>>
    BoundingBox(Args&&... args): Box(std::forward<Args>(args)...), _children(8), _contents{}, _level(0) {}
    */

    BoundingBox(const Point& lower, const Point& upper);
    explicit BoundingBox(double x0=0.0, double y0=0.0, double z0=0.0, double x1=1.0, double y1=1.0, double z1=1.0);
    BoundingBox(const BoundingBox&) = delete;
    BoundingBox(BoundingBox&&) = default;
    ~BoundingBox() = default;
    
    BoundingBox& operator=(const BoundingBox&) = delete;
    BoundingBox& operator=(BoundingBox&&) = default;

    Node& operator[](std::size_t i) noexcept{ return _children[i]; }
    const Node& operator[](std::size_t i) const noexcept{ return _children[i]; }
    
    Node& operator()(std::size_t i) noexcept{ return _contents[i]; }
    const Node& operator()(std::size_t i) const noexcept{ return _contents[i]; }
    void push(Node& node) noexcept{ _contents.push_back(std::move(node)); }
    void push(Node&& node) noexcept{ _contents.push_back(std::move(node)); }
    
    std::size_t level(){ return _level; }
    void setLevel(std::size_t level){ _level = level; }

    bool contentsOverlap(const Shape& other) const noexcept override;

    std::size_t size() const noexcept;
    constexpr std::size_t max_size() const noexcept{ return 8; }
    std::size_t numContents() const noexcept{ return _contents.size(); }
    bool empty() const noexcept;
    bool full() const noexcept;

    std::size_t octant(const Shape& other) const noexcept;

    protected:
    std::array<Node, 8> _children;
    std::vector<Node> _contents;
    std::size_t _level;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif // BOUNDING_BOX_H