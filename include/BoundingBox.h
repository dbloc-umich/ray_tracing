#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H
#include "Box.h"
#include <memory>
#include <vector>

using Node = std::unique_ptr<Shape>;
using NodeList = std::vector<Node>;

class BoundingBox: public Box{
    public:
    BoundingBox(const Point& lower, const Point& upper);
    explicit BoundingBox(double x0=0.0, double y0=0.0, double z0=0.0, double x1=1.0, double y1=1.0, double z1=1.0);
    explicit BoundingBox(NodeList& values);
    BoundingBox(const BoundingBox&) = delete;
    BoundingBox(BoundingBox&&) = default;
    ~BoundingBox() = default;
    
    BoundingBox& operator=(const BoundingBox&) = delete;
    BoundingBox& operator=(BoundingBox&&) = default;

    Node& operator[](std::size_t i) noexcept;
    const Node& operator[](std::size_t i) const noexcept{ return _children[i]; }
    std::size_t level(){ return _level; }
    void setLevel(std::size_t level){ _level = level; }

    bool contentsOverlap(const Shape& other) const noexcept override;

    std::size_t size() const noexcept{ return _children.size(); }
    constexpr std::size_t max_size() const noexcept{ return 8; }
    bool empty() const noexcept{ return size() == 0; }
    bool full() const noexcept{ return size() == 8; }

    protected:
    NodeList _children;
    std::size_t _level;
    std::ostream& print(std::ostream& os) const noexcept override;
};
#endif // BOUNDING_BOX_H