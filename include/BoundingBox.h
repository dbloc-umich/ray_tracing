#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H
#include "Box.h"
#include <memory>
#include <type_traits>
#include <vector>

using Node = std::unique_ptr<Shape>;
using NodeList = std::vector<Node>;

class BoundingBox: public Box{
    public:
    template <typename... Args,
              typename = std::enable_if_t<!std::is_base_of<BoundingBox, std::decay_t<Args>...>::value>>
    BoundingBox(Args&&... args): Box(std::forward<Args>(args)...), _children(8), _contents{}, _level(0) {}

    BoundingBox(const BoundingBox&) = delete;
    BoundingBox(BoundingBox&&) = default;
    ~BoundingBox() = default;
    
    BoundingBox& operator=(const BoundingBox&) = delete;
    BoundingBox& operator=(BoundingBox&&) = default;

    Node& operator[](std::size_t i) noexcept{ return _children[i]; }
    const Node& operator[](std::size_t i) const noexcept{ return _children[i]; }
    Node& operator()(std::size_t i) noexcept{ return _contents[i]; }
    const Node& operator()(std::size_t i) const noexcept{ return _contents[i]; }
    std::size_t level(){ return _level; }
    void setLevel(std::size_t level){ _level = level; }

    bool contentsOverlap(const Shape& other) const noexcept override;

    std::size_t size() const noexcept;
    constexpr std::size_t max_size() const noexcept{ return 8; }
    std::size_t numContents() const noexcept{ return _contents.size(); }
    bool empty() const noexcept;
    bool full() const noexcept;

    protected:
    NodeList _children, _contents;
    std::size_t _level;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif // BOUNDING_BOX_H