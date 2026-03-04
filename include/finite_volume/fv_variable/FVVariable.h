#ifndef FV_VARIABLE_H
#define FV_VARIABLE_H

#include <string>
class FVVariable{
    public:
    explicit FVVariable(const std::string& name, unsigned short dim = 1): _name(name), _dim(dim) {}
    std::string name() const noexcept{ return _name; };
    unsigned short dim() const noexcept{ return _dim; };
    virtual bool aux() const noexcept = 0;

    protected:
    std::string _name;
    unsigned short _dim;
};

#endif