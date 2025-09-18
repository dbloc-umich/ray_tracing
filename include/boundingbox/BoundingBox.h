#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H
#include "Box.h"
#include <array>
#include <memory>

using Node = std::unique_ptr<Shape>;
class BoundingBox: public Box{
    public:
    BoundingBox(const Point& lower, const Point& upper);
    explicit BoundingBox(double x0=0.0, double y0=0.0, double z0=0.0, double x1=1.0, double y1=1.0, double z1=1.0);
    BoundingBox(const BoundingBox&) = delete;
    BoundingBox(BoundingBox&&) = default;
    ~BoundingBox() = default;
    
    BoundingBox& operator=(const BoundingBox&) = delete;
    BoundingBox& operator=(BoundingBox&&) = default;

    // Access the children
    Node& operator[](std::size_t i) noexcept{ return _children[i]; }
    const Node& operator[](std::size_t i) const noexcept{ return _children[i]; }
    
    std::size_t level(){ return _level; }
    void setLevel(std::size_t level){ _level = level; }

    bool contentsOverlap(const Shape& other) const noexcept override;

    virtual double solidVolume() const noexcept;
    std::size_t size() const noexcept;
    constexpr std::size_t max_size() const noexcept{ return 8; }
    virtual bool empty() const noexcept;
    bool full() const noexcept;

    protected:
    std::array<Node, 8> _children;
    std::size_t _level;
    std::ostream& printTabs(std::ostream& os, std::size_t n) const noexcept;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif // BOUNDING_BOX_H